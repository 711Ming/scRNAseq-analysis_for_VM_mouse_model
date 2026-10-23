# -*- coding: utf-8 -*-
"""
小鼠SWS单细胞转录组分析代码

分析顺序：
1. 原始10x数据读取、质量控制、标准化、降维聚类及cluster marker识别；
2. 根据marker结果完成细胞大类注释，并提取内皮细胞及其亚群；
3. 分析内皮细胞亚群的差异表达、功能富集、细胞周期及信号通路特征；
4. 进行MELD疾病状态评分及动静脉特征评分；
5. 分析tdTomato阳性内皮细胞及其亚群的转录特征；
6. 基于loom数据进行RNA velocity分析。

输入目录：
data/matrix/  各样本10x表达矩阵目录
data/loom/    各样本loom文件

输出目录：
results/      处理后的数据、统计表及差异/富集分析结果
results/figures/  可视化结果
"""

from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import matplotlib.pyplot as plt
from scipy import sparse

# ==================== 0.参数设置 ====================
# 按样本实际存放位置修改数据目录。所有中间h5ad文件均采用通用文件名。
ROOT = Path(__file__).resolve().parents[1]
MATRIX_DIR = ROOT / "data" / "matrix"
LOOM_DIR = ROOT / "data" / "loom"
DESEQ2_DIR = ROOT / "results" / "deseq2_results"
OUTPUT_DIR = ROOT / "results"
FIGURE_DIR = OUTPUT_DIR / "figures"

# 以下为已上传Notebook中实际启用的10x矩阵样本。
# 原Notebook中WT_2和WT_3矩阵读取语句被注释，因此这里不将其列入原始矩阵样本。
SAMPLES = [
    {"sample": "CM12_1", "matrix": "CM12_1", "group": "CM12", "condition": "CM"},
    {"sample": "CM12_2", "matrix": "CM12_2", "group": "CM12", "condition": "CM"},
    {"sample": "CM12_3", "matrix": "CM12_3", "group": "CM12", "condition": "CM"},
    {"sample": "CM28_1", "matrix": "CM28_1", "group": "CM28", "condition": "CM"},
    {"sample": "CM28_2", "matrix": "CM28_2", "group": "CM28", "condition": "CM"},
    {"sample": "WT_1", "matrix": "WT28", "group": "WT", "condition": "WT"},
]

# 根据cluster marker结果填写聚类编号与细胞类型的对应关系。
# 内皮细胞亚群标签使用aEC、vEC、cEC和dEC；其余细胞类型按实际注释填写。
CLUSTER_TO_CELLTYPE = {
    # "0": "aEC",
    # "1": "vEC",
    # "2": "cEC",
    # "3": "Astrocyte",
}

EC_TYPES = ["aEC", "vEC", "cEC", "dEC"]

# 质量控制及降维聚类参数。
MIN_GENES_PER_CELL = 300
MIN_CELLS_PER_GENE = 50
MAX_GENES_PER_CELL = 4000
MAX_MT_PERCENT = 15
TARGET_SUM = 10000
N_TOP_HVG = 10000
N_NEIGHBORS = 100
N_PCS = 30
LOUVAIN_RESOLUTION = 1.0

# loom文件样本名。文件应位于data/loom/目录下。
LOOM_SAMPLES = [
    "CM12_1", "CM12_2", "CM12_3", "CM12_4",
    "CM28_1", "CM28_2", "WT_1", "WT_2", "WT_3",
]

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = str(FIGURE_DIR)
sc.settings.verbosity = 2


# 将稀疏表达矩阵中的单列基因表达转换为一维数组，便于后续状态判断和统计。
def to_dense_1d(x):
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x).reshape(-1)


# 根据当前AnnData对象中的基因名筛选可用基因，避免因数据集缺少某些基因而中断分析。
def available_genes(adata, genes, use_raw=True):
    gene_names = adata.raw.var_names if use_raw and adata.raw is not None else adata.var_names
    return [gene for gene in genes if gene in gene_names]


# 按生物学样本汇总细胞的原始整数计数，生成基因×样本的pseudo-bulk矩阵。
# 该矩阵可作为后续DESeq2分析的输入；不能用log1p后的表达量替代原始计数。
def export_pseudobulk(adata, label):
    if "counts" not in adata.layers:
        warnings.warn(f"{label}: 缺少原始计数层，跳过pseudo-bulk导出。")
        return
    if adata.n_obs == 0:
        return

    count_matrix = adata.layers["counts"]
    samples = adata.obs["sample"].astype(str).drop_duplicates().tolist()
    columns = {}
    for sample_id in samples:
        mask = adata.obs["sample"].astype(str).to_numpy() == sample_id
        summed = count_matrix[mask].sum(axis=0)
        if sparse.issparse(summed):
            summed = summed.toarray()
        columns[sample_id] = np.asarray(summed).reshape(-1)

    count_df = pd.DataFrame(columns, index=adata.var_names)
    count_df.index.name = "gene"
    count_df = count_df.round().astype("int64")
    count_df.to_csv(OUTPUT_DIR / f"{label}_pseudobulk_counts.csv")

    meta_cols = [col for col in ["sample", "group", "condition"] if col in adata.obs]
    adata.obs[meta_cols].drop_duplicates().set_index("sample").to_csv(
        OUTPUT_DIR / f"{label}_sample_info.csv"
    )


# 按原始Notebook中的做法，对单细胞评分进行组间Mann–Whitney U检验，并进行BH校正。
# 该结果反映细胞层面的分布差异；正式生物学推断还应结合样本层面的重复设计解释。
def run_cell_level_pairwise_test(adata, score_col, label):
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests

    comparisons = [("WT", "CM12"), ("WT", "CM28"), ("CM12", "CM28")]
    rows = []
    data = adata.obs[["group", score_col]].dropna()
    for group1, group2 in comparisons:
        x = data.loc[data["group"].astype(str) == group1, score_col].astype(float)
        y = data.loc[data["group"].astype(str) == group2, score_col].astype(float)
        if len(x) == 0 or len(y) == 0:
            continue
        stat, pvalue = mannwhitneyu(x, y, alternative="two-sided")
        rows.append({
            "score": score_col,
            "comparison": f"{group1}_vs_{group2}",
            "n_group1_cells": len(x),
            "n_group2_cells": len(y),
            "mean_group1": x.mean(),
            "mean_group2": y.mean(),
            "median_group1": x.median(),
            "median_group2": y.median(),
            "U_statistic": stat,
            "p_value": pvalue,
        })
    if rows:
        result = pd.DataFrame(rows)
        result["p_adj_BH"] = multipletests(result["p_value"], method="fdr_bh")[1]
        result.to_csv(OUTPUT_DIR / f"{label}_cell_level_tests.csv", index=False)


# 对给定差异基因列表进行GO Biological Process和KEGG过度富集分析。
def run_enrichment(gene_list, prefix):
    import gseapy as gp

    gene_list = list(pd.Series(gene_list).dropna().astype(str).unique())
    if not gene_list:
        return
    for gene_set, suffix in [
        ("GO_Biological_Process_2021", "GO_BP"),
        ("KEGG_2016", "KEGG"),
    ]:
        result = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_set,
            organism="Mouse",
            cutoff=0.05,
            outdir=None,
        )
        result.res2d.to_csv(OUTPUT_DIR / f"{prefix}_{suffix}_enrichment.csv", index=False)


# ==================== 1.原始数据处理、质量控制及聚类 ====================
# 逐个读取10x矩阵，添加样本、组别和条件信息，并按共同基因集合合并样本。
# 随后过滤低质量细胞和低表达基因，保存原始counts，再进行标准化和降维聚类。
# 本函数同时输出QC统计表、UMAP、候选marker图及各cluster的marker基因表。
def load_and_cluster():
    # 1.1 读取10x矩阵并为细胞条形码添加样本后缀，防止不同样本的barcode重复。
    sample_objects = []
    for item in SAMPLES:
        matrix_path = MATRIX_DIR / item["matrix"]
        sample = sc.read_10x_mtx(
            str(matrix_path),
            var_names="gene_symbols",
            make_unique=True,
        )
        sample.var_names_make_unique()
        sample.obs["sample"] = item["sample"]
        sample.obs["group"] = item["group"]
        sample.obs["condition"] = item["condition"]
        sample.obs_names = [f"{barcode}{item['sample']}" for barcode in sample.obs_names.astype(str)]
        sample_objects.append(sample)

    # 1.2 仅保留所有样本共有的基因，合并成统一AnnData对象。
    adata = ad.concat(sample_objects, join="inner", index_unique=None, merge="same")
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    # 原始Notebook中使用线粒体比例及检测基因数进行过滤；上传文件未见Scrublet调用。
    # 若实际分析使用了Scrublet，应将真实运行代码补入此处并确保论文方法与代码一致。
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # 1.3 汇总每个样本的细胞数、检测基因数、UMI数及线粒体比例。
    qc = (
        adata.obs.groupby("sample", observed=True)
        .agg(
            n_cells=("sample", "size"),
            median_genes=("n_genes_by_counts", "median"),
            median_UMI=("total_counts", "median"),
            median_mt_percent=("pct_counts_mt", "median"),
        )
        .reset_index()
    )
    qc.to_csv(OUTPUT_DIR / "QC_summary_before_final_filter.csv", index=False)
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save="_QC.pdf",
    )

    keep = (
        (adata.obs["pct_counts_mt"] < MAX_MT_PERCENT)
        & (adata.obs["n_genes_by_counts"] < MAX_GENES_PER_CELL)
    )
    adata = adata[keep].copy()

    # 1.4 在标准化前保存原始整数计数，供pseudo-bulk计数矩阵构建使用。
    adata.layers["counts"] = adata.X.copy()
    # 1.5 每细胞总计数归一化至10,000并进行log1p转换。
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    # 1.6 筛选高变基因，并在高变基因空间进行PCA、Harmony批次校正及无监督聚类。
    sc.pp.highly_variable_genes(
        adata,
        flavor="cell_ranger",
        n_top_genes=min(N_TOP_HVG, adata.n_vars),
    )
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.tl.pca(adata_hvg, svd_solver="arpack")
    sce.pp.harmony_integrate(adata_hvg, key="sample")

    n_neighbors = min(N_NEIGHBORS, max(2, adata_hvg.n_obs - 1))
    n_pcs = min(N_PCS, adata_hvg.obsm["X_pca_harmony"].shape[1])
    sc.pp.neighbors(
        adata_hvg,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep="X_pca_harmony",
    )
    sc.tl.umap(adata_hvg)
    sc.tl.louvain(
        adata_hvg,
        resolution=LOUVAIN_RESOLUTION,
        key_added="louvain_1.0",
    )

    adata.obs["louvain_1.0"] = adata_hvg.obs["louvain_1.0"].reindex(
        adata.obs_names
    ).astype("category")
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"].copy()
    adata.obsm["X_pca_harmony"] = adata_hvg.obsm["X_pca_harmony"].copy()

    sc.tl.rank_genes_groups(
        adata,
        groupby="louvain_1.0",
        method="wilcoxon",
        tie_correct=True,
        use_raw=True,
    )
    sc.get.rank_genes_groups_df(adata, group=None).to_csv(
        OUTPUT_DIR / "cluster_marker_genes.csv",
        index=False,
    )

    marker_genes = [
        "Pecam1", "Cdh5", "Vwf", "Acta2", "Tagln", "Myh11",
        "Pdgfrb", "Rgs5", "Cspg4", "Des", "Kcnj8", "Abcc9",
        "Col1a1", "Col1a2", "Dcn", "Gfap", "Aqp4", "Pdgfra",
        "Plp1", "Mbp", "Tmem119", "Gpr34", "Ptprc",
    ]
    marker_genes = available_genes(adata, marker_genes)
    if marker_genes:
        sc.pl.umap(
            adata,
            color=marker_genes,
            ncols=4,
            frameon=False,
            save="_candidate_markers.pdf",
        )
        sc.pl.dotplot(
            adata,
            marker_genes,
            groupby="louvain_1.0",
            use_raw=True,
            standard_scale="var",
            save="_cluster_marker_dotplot.pdf",
        )

    adata.write_h5ad(OUTPUT_DIR / "all_cells_clustered.h5ad", compression="gzip")
    return adata


# ==================== 2.细胞大类注释 ====================
# 根据cluster marker表和候选marker可视化结果，人工确认cluster与细胞类型的对应关系。
# 由于原Notebook没有包含完整的cluster-to-celltype映射，CLUSTER_TO_CELLTYPE需根据实际注释填写。
def annotate_cells(adata):
    clusters = set(adata.obs["louvain_1.0"].astype(str).unique())
    mapping = {str(key): value for key, value in CLUSTER_TO_CELLTYPE.items()}
    missing = sorted(clusters - set(mapping))
    if missing:
        raise SystemExit(
            "请根据cluster_marker_genes.csv和marker图填写CLUSTER_TO_CELLTYPE；"
            f"尚未注释的cluster：{missing}。聚类对象已保存至results/all_cells_clustered.h5ad。"
        )

    adata.obs["cell_type_large"] = (
        adata.obs["louvain_1.0"].astype(str).map(mapping).astype("category")
    )
    adata.write_h5ad(OUTPUT_DIR / "all_cells_annotated.h5ad", compression="gzip")

    sc.pl.umap(
        adata,
        color=["cell_type_large", "group"],
        wspace=0.35,
        save="_cell_type_large.pdf",
    )

    counts = (
        adata.obs.groupby(["sample", "group", "condition", "cell_type_large"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    counts.to_csv(OUTPUT_DIR / "celltype_counts_by_sample.csv", index=False)
    pd.crosstab(adata.obs["group"], adata.obs["cell_type_large"]).to_csv(
        OUTPUT_DIR / "celltype_counts_by_group.csv"
    )
    return adata


# ==================== 3.tdTomato状态判定 ====================
# 依据标准化前保存的tdTomato原始计数判断阳性状态；表达计数大于0判定为阳性。
# 如果表达矩阵中不存在tdTomato基因，则将状态标记为Unknown，不将未匹配细胞归为阴性。
def assign_tdTomato_status(adata):
    if "tdTomato" not in adata.var_names or "counts" not in adata.layers:
        adata.obs["tdTomato_status"] = "Unknown"
        adata.obs["tdTomato_pos"] = False
        return adata

    gene_index = adata.var_names.get_loc("tdTomato")
    counts = adata.layers["counts"][:, gene_index]
    expression = to_dense_1d(counts)
    status = np.where(expression > 0, "Positive", "Negative")
    adata.obs["tdTomato_status"] = pd.Categorical(
        status,
        categories=["Negative", "Positive", "Unknown"],
    )
    adata.obs["tdTomato_pos"] = adata.obs["tdTomato_status"].astype(str).eq("Positive")
    pd.crosstab(adata.obs["sample"], adata.obs["tdTomato_status"]).to_csv(
        OUTPUT_DIR / "tdTomato_by_sample.csv"
    )
    pd.crosstab(adata.obs["cell_type_large"], adata.obs["tdTomato_status"]).to_csv(
        OUTPUT_DIR / "tdTomato_by_celltype.csv"
    )
    sc.pl.umap(adata, color=["tdTomato_status", "cell_type_large"], save="_tdTomato.pdf")
    return adata


# ==================== 4.去除fEC并提取内皮细胞 ====================
# 按已注释的cell_type_large去除fEC，保留其余目标血管相关细胞；
# 随后提取aEC、vEC、cEC和dEC，输出样本及亚群细胞数量统计。
def extract_cell_subsets(adata):
    adata = assign_tdTomato_status(adata)
    vascular = adata[~adata.obs["cell_type_large"].astype(str).isin(["fEC"])].copy()
    vascular.write_h5ad(OUTPUT_DIR / "vascular_cells_no_fEC.h5ad", compression="gzip")

    ec = vascular[vascular.obs["cell_type_large"].astype(str).isin(EC_TYPES)].copy()
    if ec.n_obs == 0:
        raise ValueError(f"未提取到EC细胞，请检查注释标签是否为{EC_TYPES}。")
    ec.write_h5ad(OUTPUT_DIR / "endothelial_cells.h5ad", compression="gzip")

    ec_counts = (
        ec.obs.groupby(["sample", "group", "condition", "cell_type_large"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    ec_counts.to_csv(OUTPUT_DIR / "EC_subtype_counts_by_sample.csv", index=False)
    return vascular, ec


# ==================== 5.内皮细胞亚群详细分析 ====================
# 5.1 检查内皮细胞身份标志及动静脉、血管生成和信号通路相关基因表达；
# 5.2 筛选dEC相对于其他内皮细胞的marker，并对上调marker进行功能富集；
# 5.3 根据S期和G2/M期基因集计算细胞周期评分，统计增殖细胞比例；
# 5.4 计算vEC静脉特征评分和aEC动脉特征评分；
# 5.5 导出各内皮亚群的pseudo-bulk原始计数矩阵，供DESeq2分析。
def analyze_endothelial_cells(ec):
    ec_markers = [
        "Pecam1", "Cdh5", "Vwf", "Slco1c1", "Cldn5", "Tjp1", "Mfsd2a",
        "Nr2f2", "Ephb4", "Nrp2", "Ackr1", "Emcn", "Bmp4", "Apln",
        "Efnb2", "Dll4", "Notch1", "Hey1", "Hey2", "Gja5", "Sox17", "Alk1",
        "Vegfa", "Vegfb", "Vegfc", "Kdr", "Flt1", "Angpt1", "Angpt2", "Tek", "Tie1", "Foxo1", "Ptprb",
        "Plcd1", "Plcb1", "Prkca", "Pik3r3", "Pik3ca", "Pik3r1", "Pdk1", "Akt1", "Rptor", "Mtor", "Eif4e", "Eif4ebp1",
        "Egfr", "Grb2", "Sos1", "Kras", "Braf", "Map2k1", "Mapk1", "Mapk3", "Fos", "Jun", "Egr1",
        "Mki67", "Top2a", "Mcm2", "Cdk1", "Birc5",
    ]
    ec_markers = available_genes(ec, ec_markers)
    if ec_markers:
        sc.pl.dotplot(
            ec,
            ec_markers,
            groupby=["cell_type_large", "group"],
            use_raw=True,
            standard_scale="var",
            save="_EC_marker_pathway_dotplot.pdf",
        )
        sc.pl.umap(
            ec,
            color=["Cdh5", "Slco1c1", "cell_type_large"],
            save="_EC_identity.pdf",
        )

    if "dEC" in ec.obs["cell_type_large"].astype(str).unique():
        sc.tl.rank_genes_groups(
            ec,
            groupby="cell_type_large",
            groups=["dEC"],
            reference="rest",
            method="wilcoxon",
            tie_correct=True,
            use_raw=True,
        )
        dec = sc.get.rank_genes_groups_df(ec, group="dEC")
        dec.to_csv(OUTPUT_DIR / "dEC_markers_all.csv", index=False)
        dec_up = dec[
            (dec["pvals_adj"] < 0.05)
            & (dec["logfoldchanges"] > 0.25)
        ].copy()
        dec_up.to_csv(OUTPUT_DIR / "dEC_markers_up.csv", index=False)
        try:
            run_enrichment(dec_up["names"].tolist(), "dEC_up")
        except Exception as exc:
            warnings.warn(f"dEC富集分析未完成：{exc}")

    s_genes = [
        "Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung",
        "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Hells", "Rfc2",
        "Rpa2", "Nasp", "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7",
        "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1",
        "Tipin", "Dscc1", "Blm", "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b",
        "Brip1", "E2f8",
    ]
    g2m_genes = [
        "Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80",
        "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", "Tacc3", "Smc4",
        "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", "Bub1", "Kif11", "Anln", "Lbr",
        "Cenpe", "Ctcf", "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa",
    ]
    s_genes = available_genes(ec, s_genes, use_raw=False)
    g2m_genes = available_genes(ec, g2m_genes, use_raw=False)
    if s_genes and g2m_genes:
        sc.tl.score_genes(ec, gene_list=s_genes, score_name="S_phase_score", use_raw=False)
        sc.tl.score_genes(ec, gene_list=g2m_genes, score_name="G2M_phase_score", use_raw=False)
        ec.obs["CellCyclePhase"] = [
            "Quiescent" if (s + g2m) < 0.2 else "Proliferating"
            for s, g2m in zip(ec.obs["S_phase_score"], ec.obs["G2M_phase_score"])
        ]
        ec.obs["CellCyclePhase"] = pd.Categorical(
            ec.obs["CellCyclePhase"],
            categories=["Quiescent", "Proliferating"],
        )
        phase_summary = (
            ec.obs.groupby(["group", "cell_type_large", "CellCyclePhase"], observed=True)
            .size()
            .unstack(fill_value=0)
            .reset_index()
        )
        for column in ["Quiescent", "Proliferating"]:
            if column not in phase_summary:
                phase_summary[column] = 0
        denominator = phase_summary["Quiescent"] + phase_summary["Proliferating"]
        phase_summary["Proliferating_ratio"] = np.divide(
            phase_summary["Proliferating"],
            denominator,
            out=np.full(len(phase_summary), np.nan),
            where=denominator.to_numpy() > 0,
        )
        phase_summary.to_csv(OUTPUT_DIR / "EC_proliferation_summary.csv", index=False)
        sc.pl.umap(ec, color=["CellCyclePhase", "group"], save="_EC_cell_cycle.pdf")

    venous_genes = ["Nr2f2", "Ephb4", "Nrp2", "Ackr1", "Emcn", "Bmp4", "Alk3", "Apln"]
    arterial_genes = ["Efnb2", "Dll4", "Notch1", "Hey1", "Hey2", "Gja5", "Sox17", "Alk1"]
    if "vEC" in ec.obs["cell_type_large"].astype(str).unique():
        vec = ec[ec.obs["cell_type_large"].astype(str) == "vEC"].copy()
        genes = available_genes(vec, venous_genes, use_raw=False)
        if genes:
            sc.tl.score_genes(vec, gene_list=genes, score_name="venous_score", ctrl_size=50, use_raw=False)
            ec.obs.loc[vec.obs_names, "venous_score"] = vec.obs["venous_score"].to_numpy()
            run_cell_level_pairwise_test(vec, "venous_score", "vEC_venous_score")
    if "aEC" in ec.obs["cell_type_large"].astype(str).unique():
        aec = ec[ec.obs["cell_type_large"].astype(str) == "aEC"].copy()
        genes = available_genes(aec, arterial_genes, use_raw=False)
        if genes:
            sc.tl.score_genes(aec, gene_list=genes, score_name="arterial_score", ctrl_size=50, use_raw=False)
            ec.obs.loc[aec.obs_names, "arterial_score"] = aec.obs["arterial_score"].to_numpy()
            run_cell_level_pairwise_test(aec, "arterial_score", "aEC_arterial_score")

    score_columns = [col for col in ["venous_score", "arterial_score"] if col in ec.obs]
    if score_columns:
        ec.obs[["sample", "group", "condition", "cell_type_large"] + score_columns].groupby(
            ["sample", "group", "condition", "cell_type_large"], observed=True
        ).mean(numeric_only=True).reset_index().to_csv(
            OUTPUT_DIR / "EC_arterial_venous_scores_by_sample.csv", index=False
        )

    for cell_type in EC_TYPES:
        export_pseudobulk(
            ec[ec.obs["cell_type_large"].astype(str) == cell_type].copy(),
            f"{cell_type}_EC",
        )

    ec.write_h5ad(OUTPUT_DIR / "endothelial_cells_analyzed.h5ad", compression="gzip")
    return ec


# ==================== 6.信号通路及血管生成相关基因分析 ====================
# 按原Notebook展示PI3K–AKT–mTOR、MAPK、VEGF及ANGPT–TIE相关基因，
# 并统计Angpt1、Angpt2、Tek、Tie1和Foxo1在不同内皮亚群/组别中的平均表达及阳性比例。
def analyze_pathway_and_angiogenesis(ec):
    pathway_sets = {
        "PI3K_AKT_mTOR": [
            "Plcd1", "Plcb1", "Prkca", "Pik3r3", "Pik3ca", "Pik3r1",
            "Pdk1", "Akt1", "Rptor", "Mtor", "Eif4e", "Eif4ebp1",
            "Vegfa", "Angpt2",
        ],
        "MAPK": [
            "Egfr", "Grb2", "Sos1", "Kras", "Braf", "Map2k1",
            "Mapk1", "Mapk3", "Fos", "Jun", "Egr1",
        ],
        "VEGF": ["Vegfa", "Vegfb", "Vegfc", "Kdr", "Flt1", "Nrp1", "Nrp2"],
        "ANGPT_TIE": ["Angpt1", "Angpt2", "Tek", "Tie1", "Foxo1", "Ptprb"],
    }

    for set_name, genes in pathway_sets.items():
        genes_use = available_genes(ec, genes, use_raw=True)
        if not genes_use:
            continue
        sc.pl.dotplot(
            ec,
            genes_use,
            groupby=["cell_type_large", "group"],
            use_raw=True,
            standard_scale="var",
            save=f"_{set_name}_dotplot.pdf",
        )

    target_genes = ["Angpt1", "Angpt2", "Tek", "Tie1", "Foxo1", "Ptprb"]
    target_genes = available_genes(ec, target_genes, use_raw=False)
    if not target_genes:
        return

    # X为log1p标准化表达矩阵；表达值大于0用于统计检测到该基因的细胞比例。
    expr = ec[:, target_genes].X
    if sparse.issparse(expr):
        expr = expr.toarray()
    expr_df = pd.DataFrame(expr, index=ec.obs_names, columns=target_genes)
    expr_df["sample"] = ec.obs["sample"].astype(str).values
    expr_df["group"] = ec.obs["group"].astype(str).values
    expr_df["cell_type_large"] = ec.obs["cell_type_large"].astype(str).values

    summary_rows = []
    for (cell_type, group), sub in expr_df.groupby(
        ["cell_type_large", "group"], observed=True
    ):
        for gene in target_genes:
            values = sub[gene].astype(float)
            summary_rows.append({
                "cell_type_large": cell_type,
                "group": group,
                "gene": gene,
                "mean_expression": values.mean(),
                "positive_cell_ratio": (values > 0).mean(),
                "n_cells": len(values),
            })
    pd.DataFrame(summary_rows).to_csv(
        OUTPUT_DIR / "EC_angiogenic_gene_expression_summary.csv",
        index=False,
    )

    # 导出逐细胞表达表，便于后续绘图或按样本重新汇总。
    expr_df.to_csv(OUTPUT_DIR / "EC_angiogenic_gene_expression_by_cell.csv")

    # 对原Notebook重点关注的Lrg1、Tek、Vegfc等基因，比较aEC和vEC在不同组别中的表达。
    response_genes = [
        "Lrg1", "Tek", "Vegfc", "Klf2", "Fn1",
        "Igfbp3", "Egfl7", "Mgp", "Mmrn2", "Angpt2",
    ]
    response_genes = available_genes(ec, response_genes, use_raw=False)
    ec_sub = ec[
        ec.obs["cell_type_large"].astype(str).isin(["aEC", "vEC"])
    ].copy()
    if response_genes and ec_sub.n_obs > 0:
        ec_sub.obs["celltype_condition"] = (
            ec_sub.obs["cell_type_large"].astype(str)
            + "_"
            + ec_sub.obs["condition"].astype(str)
        )
        sc.pl.violin(
            ec_sub,
            keys=response_genes,
            groupby="celltype_condition",
            rotation=90,
            stripplot=False,
            save="_aEC_vEC_response_genes.pdf",
        )
        response_expr = ec_sub[:, response_genes].X
        if sparse.issparse(response_expr):
            response_expr = response_expr.toarray()
        response_df = pd.DataFrame(
            response_expr,
            index=ec_sub.obs_names,
            columns=response_genes,
        )
        response_df["sample"] = ec_sub.obs["sample"].astype(str).values
        response_df["condition"] = ec_sub.obs["condition"].astype(str).values
        response_df["cell_type_large"] = ec_sub.obs["cell_type_large"].astype(str).values
        response_df.to_csv(
            OUTPUT_DIR / "aEC_vEC_response_gene_expression_by_cell.csv"
        )





# ==================== 7.tdTomato阳性内皮细胞分析 ====================
# 统计tdTomato阳性细胞在样本、实验组及内皮亚群中的分布；
# 对阳性内皮细胞及aEC、vEC子集进行组间差异表达分析，并比较两类亚群的响应差异。
def analyze_tdTomato(ec):
    if "tdTomato_status" not in ec.obs:
        return
    positive = ec[ec.obs["tdTomato_status"].astype(str) == "Positive"].copy()
    if positive.n_obs == 0:
        return

    pd.crosstab(positive.obs["condition"], positive.obs["cell_type_large"]).to_csv(
        OUTPUT_DIR / "tdTomato_positive_EC_by_condition.csv"
    )
    export_pseudobulk(positive, "tdTomato_positive_EC")

    if positive.obs["condition"].astype(str).nunique() >= 2:
        sc.tl.rank_genes_groups(
            positive,
            groupby="condition",
            groups=["CM"],
            reference="WT",
            method="wilcoxon",
            tie_correct=True,
            use_raw=True,
        )
        deg = sc.get.rank_genes_groups_df(positive, group="CM")
        deg.to_csv(OUTPUT_DIR / "tdTomato_positive_EC_CM_vs_WT_DEG.csv", index=False)

        try:
            import gseapy as gp
            ranked = deg[["names", "scores"]].dropna().drop_duplicates("names")
            if len(ranked):
                ranked.columns = ["gene", "score"]
                prerank = gp.prerank(
                    rnk=ranked,
                    gene_sets="GO_Biological_Process_2021",
                    permutation_num=1000,
                    min_size=10,
                    max_size=300,
                    outdir=None,
                )
                prerank.res2d.to_csv(OUTPUT_DIR / "tdTomato_positive_EC_CM_vs_WT_GSEA.csv", index=False)
            sig = deg[deg["pvals_adj"] < 0.05]
            run_enrichment(sig.loc[sig["logfoldchanges"] > 0, "names"], "tdTomato_positive_EC_up")
            run_enrichment(sig.loc[sig["logfoldchanges"] < 0, "names"], "tdTomato_positive_EC_down")
        except Exception as exc:
            warnings.warn(f"tdTomato阳性EC富集分析未完成：{exc}")

    subtype_results = {}
    for cell_type in ["aEC", "vEC"]:
        subset = positive[positive.obs["cell_type_large"].astype(str) == cell_type].copy()
        if subset.n_obs == 0 or subset.obs["condition"].astype(str).nunique() < 2:
            continue
        sc.tl.rank_genes_groups(
            subset,
            groupby="condition",
            groups=["CM"],
            reference="WT",
            method="wilcoxon",
            tie_correct=True,
            use_raw=True,
        )
        result = sc.get.rank_genes_groups_df(subset, group="CM")
        result.to_csv(OUTPUT_DIR / f"tdTomato_positive_{cell_type}_CM_vs_WT_DEG.csv", index=False)
        subtype_results[cell_type] = result

    if "aEC" in subtype_results and "vEC" in subtype_results:
        a = subtype_results["aEC"][["names", "logfoldchanges", "pvals_adj"]].rename(
            columns={"logfoldchanges": "logFC_aEC", "pvals_adj": "padj_aEC"}
        )
        v = subtype_results["vEC"][["names", "logfoldchanges", "pvals_adj"]].rename(
            columns={"logfoldchanges": "logFC_vEC", "pvals_adj": "padj_vEC"}
        )
        comparison = a.merge(v, on="names", how="outer")
        comparison["delta_logFC_vEC_minus_aEC"] = comparison["logFC_vEC"] - comparison["logFC_aEC"]
        comparison.sort_values("delta_logFC_vEC_minus_aEC", ascending=False).to_csv(
            OUTPUT_DIR / "tdTomato_positive_aEC_vEC_response_comparison.csv",
            index=False,
        )

        both = positive[positive.obs["cell_type_large"].astype(str).isin(["aEC", "vEC"])].copy()
        def mean_expression(mask):
            x = both[mask].X
            if sparse.issparse(x):
                return np.asarray(x.mean(axis=0)).reshape(-1)
            return np.asarray(x).mean(axis=0)
        masks = {
            "CM_vEC": (both.obs["condition"].astype(str) == "CM") & (both.obs["cell_type_large"].astype(str) == "vEC"),
            "WT_vEC": (both.obs["condition"].astype(str) == "WT") & (both.obs["cell_type_large"].astype(str) == "vEC"),
            "CM_aEC": (both.obs["condition"].astype(str) == "CM") & (both.obs["cell_type_large"].astype(str) == "aEC"),
            "WT_aEC": (both.obs["condition"].astype(str) == "WT") & (both.obs["cell_type_large"].astype(str) == "aEC"),
        }
        if all(mask.sum() > 0 for mask in masks.values()):
            logfc_v = mean_expression(masks["CM_vEC"]) - mean_expression(masks["WT_vEC"])
            logfc_a = mean_expression(masks["CM_aEC"]) - mean_expression(masks["WT_aEC"])
            delta = pd.DataFrame({
                "gene": both.var_names,
                "logFC_vEC_CM_minus_WT": logfc_v,
                "logFC_aEC_CM_minus_WT": logfc_a,
                "delta_vEC_minus_aEC": logfc_v - logfc_a,
            }).sort_values("delta_vEC_minus_aEC", ascending=False)
            delta.to_csv(OUTPUT_DIR / "tdTomato_positive_aEC_vEC_delta_response.csv", index=False)
            try:
                import gseapy as gp
                ranking = delta[["gene", "delta_vEC_minus_aEC"]].dropna()
                prerank = gp.prerank(
                    rnk=ranking,
                    gene_sets="GO_Biological_Process_2021",
                    permutation_num=1000,
                    min_size=10,
                    max_size=300,
                    outdir=None,
                )
                prerank.res2d.to_csv(OUTPUT_DIR / "tdTomato_positive_aEC_vEC_delta_GSEA.csv", index=False)
            except Exception as exc:
                warnings.warn(f"aEC/vEC差异响应GSEA未完成：{exc}")




# ==================== 8.MELD疾病状态评分 ====================
# 排除dEC后，分别比较CM12与WT、CM28与CM12。
# 使用GO内皮相关基因及炎症、内皮功能障碍和血管生成相关基因构建特征集合，
# 计算细胞的组别似然评分，并按内皮亚群汇总。

def run_meld(ec):
    try:
        import meld
        import gseapy as gp
    except ImportError:
        return

    # MELD使用的基因集合包括炎症、内皮功能障碍和血管生成相关基因。
    inflammation_genes = [
        "Icam1", "Vcam1", "Sele", "Selp", "Ccl2", "Cxcl1", "Cxcl8", "Cxcl10", "Cxcl12", "Il1b",
        "Il6", "Tnf", "Tnfaip3", "Nfkb1", "Rela", "Relb", "Stat1", "Stat3", "Irf1", "Nos2",
        "Ptgs2", "Tlr2", "Tlr4", "Tlr6", "Tlr7", "Traf1", "Traf6", "Birc3", "Casp1", "Casp4",
        "Casp8", "Ccl5", "Ccr1", "Ccr5", "Cxcr4", "Cxcr2", "Cd14", "Cd40", "Cd44", "Cd86",
        "Myd88", "Nlrp3", "Socs3", "Gbp2", "Ifngr1", "Ifngr2",
    ]
    dysfunction_genes = [
        "Nos3", "Edn1", "Vcam1", "Icam1", "Sele", "Thbd", "Ptgs2", "Hmox1", "Cav1", "Klf2",
        "Klf4", "Angpt2", "Tek", "Cdh5", "Vwf", "Fn1", "Cldn5", "Ocln", "F11r", "Aqp1",
        "Nox4", "Acvrl1", "Eng", "Tie1", "Itga5", "Itgb3", "Itgav", "Rhoj", "Rock2", "Mapk14",
        "Pecam1", "Vcl", "Src", "Akt1", "Sirt1", "Dpp4", "Pdgfb", "Ece1", "Smad6", "Tgfbr2",
    ]
    angiogenesis_genes = [
        "Vegfa", "Vegfb", "Vegfc", "Kdr", "Flt1", "Angpt1", "Angpt2", "Tek", "Dll4", "Notch1",
        "Fgf2", "Fgfr1", "Fgfr2", "Fgfr3", "Cd34", "Pecam1", "Cdh5", "Itgb1", "Nrp1", "Nrp2",
        "Hif1a", "Hif3a", "Ptk2", "Robo4", "Efnb2", "Ephb4", "Tie1", "Mmp2", "Mmp9", "Mmp14",
        "Mmp19", "S1pr1", "Cxcr4", "Slit2", "Plxnd1", "Apln", "Apj", "Il8", "Bmp4", "Bmp2",
        "Timp2", "Timp3", "Col4a1", "Col18a1", "Esm1", "Adm", "Tgfb1", "Sox17", "Ets1", "Spred1",
    ]

    try:
        go_library = gp.get_library(name="GO_Biological_Process_2023", organism="Mouse")
        endothelial_genes = set()
        for term, genes in go_library.items():
            if "endothelial" in term.lower():
                endothelial_genes.update(str(g).capitalize() for g in genes)
        combined = endothelial_genes | set(dysfunction_genes + inflammation_genes + angiogenesis_genes)
    except Exception:
        combined = set(dysfunction_genes + inflammation_genes + angiogenesis_genes)

    ec_filter = ec[~ec.obs["cell_type_large"].astype(str).isin(["dEC"])].copy()
    comparisons = [("CM12", "WT", "CM12_vs_WT"), ("CM28", "CM12", "CM28_vs_CM12")]
    for case, control, label in comparisons:
        subset = ec_filter[ec_filter.obs["group"].astype(str).isin([case, control])].copy()
        if subset.n_obs == 0 or subset.obs["group"].astype(str).nunique() < 2:
            continue
        genes = [gene for gene in combined if gene in subset.var_names]
        if len(genes) < 2:
            continue
        density = meld.MELD().fit_transform(subset[:, genes].X, subset.obs["group"].astype(str))
        likelihood = meld.utils.normalize_densities(density)
        case_score = likelihood[case]
        subset.obs[f"{case}_MELD"] = np.asarray(case_score)
        summary = (
            subset.obs.groupby(["cell_type_large", "group"], observed=True)[f"{case}_MELD"]
            .mean().reset_index()
        )
        summary.to_csv(OUTPUT_DIR / f"MELD_{label}_by_celltype.csv", index=False)




# ==================== 9.读取DESeq2结果并进行功能富集 ====================
# 上传的Notebook包含pseudo-bulk矩阵导出和DESeq2结果读取/富集代码，
# 但未包含DESeq2的R执行脚本。此处读取独立DESeq2分析生成的结果表。
# 结果表需包含gene_name、log2FoldChange、pvalue、baseMean等列。
def analyze_deseq2_results():
    if not DESEQ2_DIR.exists():
        warnings.warn(
            f"未找到DESeq2结果目录：{DESEQ2_DIR}。"
            "已导出的pseudo-bulk矩阵可在DESeq2分析后放入该目录。"
        )
        return

    result_files = sorted(DESEQ2_DIR.glob("*DESeq2*.csv"))
    if not result_files:
        warnings.warn(f"{DESEQ2_DIR}中没有匹配*DESeq2*.csv的结果文件。")
        return

    import gseapy as gp

    for file_path in result_files:
        deg = pd.read_csv(file_path)
        gene_col = "gene_name" if "gene_name" in deg.columns else (
            "gene" if "gene" in deg.columns else None
        )
        required = {"log2FoldChange", "pvalue", "baseMean"}
        if gene_col is None or not required.issubset(deg.columns):
            warnings.warn(f"{file_path.name}缺少所需列，跳过。")
            continue

        deg = deg.dropna(subset=[gene_col, "log2FoldChange", "pvalue"]).copy()
        deg[gene_col] = deg[gene_col].astype(str)
        sig = deg[(deg["pvalue"] < 0.05) & (deg["baseMean"] >= 5)].copy()
        up = sig.loc[sig["log2FoldChange"] > 0, gene_col].dropna().unique().tolist()
        down = sig.loc[sig["log2FoldChange"] < 0, gene_col].dropna().unique().tolist()

        for gene_list, direction in [(up, "up"), (down, "down")]:
            if not gene_list:
                continue
            for library, library_name in [
                ("GO_Biological_Process_2021", "GO_BP"),
                ("KEGG_2016", "KEGG"),
            ]:
                try:
                    enrichment = gp.enrichr(
                        gene_list=gene_list,
                        gene_sets=library,
                        organism="Mouse",
                        cutoff=0.05,
                        outdir=None,
                    )
                    enrichment.res2d.to_csv(
                        OUTPUT_DIR / f"{file_path.stem}_{direction}_{library_name}_enrichment.csv",
                        index=False,
                    )
                except Exception as exc:
                    warnings.warn(
                        f"{file_path.name}的{direction}基因{library_name}富集未完成：{exc}"
                    )

        # GSEA按log2FoldChange排序，保留全体可用基因，不仅限于显著差异基因。
        ranked = deg[[gene_col, "log2FoldChange"]].dropna().drop_duplicates(gene_col)
        ranked = ranked.sort_values("log2FoldChange", ascending=False)
        if len(ranked) >= 10:
            ranked.columns = ["gene", "score"]
            try:
                prerank = gp.prerank(
                    rnk=ranked,
                    gene_sets="GO_Biological_Process_2021",
                    permutation_num=1000,
                    min_size=10,
                    max_size=300,
                    seed=42,
                    outdir=None,
                    verbose=False,
                )
                prerank.res2d.to_csv(
                    OUTPUT_DIR / f"{file_path.stem}_GSEA.csv",
                    index=False,
                )
            except Exception as exc:
                warnings.warn(f"{file_path.name}的GSEA未完成：{exc}")




# ==================== 10.RNA velocity分析 ====================
# 读取各样本loom文件，合并spliced、unspliced及ambiguous转录本信息，
# 与内皮细胞对象按细胞和基因交集匹配，再计算RNA velocity、velocity graph并绘制流线图。
# 原Notebook实际执行了velocity及velocity_graph；未见recover_dynamics和latent_time的执行代码。
def run_velocity(ec):
    try:
        import scvelo as scv
    except ImportError:
        return

    # 8.1 读取loom文件并添加样本信息。
    loom_adatas = []
    for sample_id in LOOM_SAMPLES:
        loom_path = LOOM_DIR / f"{sample_id}.loom"
        if not loom_path.exists():
            continue
        loom = scv.read(str(loom_path), cache=True)
        loom.obs["sample"] = sample_id
        loom.obs_names = [f"{cell}{sample_id}" for cell in loom.obs_names.astype(str)]
        loom.var_names_make_unique()
        loom_adatas.append(loom)
    if not loom_adatas:
        return

    # 8.2 合并各样本loom对象，再与EC对象按细胞和基因交集匹配。
    merged_loom = ad.concat(
        loom_adatas,
        join="outer",
        label="loom_batch",
        fill_value=0,
    )
    merged_loom.obs_names = [
        name.split(":")[-1].replace("x", "").replace("X", "")
        for name in merged_loom.obs_names.astype(str)
    ]
    merged_loom.obs_names_make_unique()
    ec.obs_names_make_unique()
    common_cells = ec.obs_names.intersection(merged_loom.obs_names)
    if len(common_cells) == 0:
        warnings.warn("loom与EC对象没有匹配到共同细胞，跳过RNA velocity。")
        return
    common_genes = ec.var_names.intersection(merged_loom.var_names)
    if len(common_genes) == 0:
        warnings.warn("loom与EC对象没有匹配到共同基因，跳过RNA velocity。")
        return

    ec_velocity = ec[common_cells, common_genes].copy()
    loom_aligned = merged_loom[ec_velocity.obs_names, common_genes].copy()
    for layer in ["spliced", "unspliced", "ambiguous"]:
        if layer in loom_aligned.layers:
            ec_velocity.layers[layer] = loom_aligned.layers[layer].copy()

    # 8.3 对匹配后的spliced/unspliced数据进行过滤、归一化及邻域矩计算。
    scv.pp.filter_and_normalize(ec_velocity)
    scv.pp.moments(ec_velocity)
    # 8.4 计算RNA velocity和细胞状态转移图，并投射至UMAP空间。
    scv.tl.velocity(ec_velocity)
    scv.tl.velocity_graph(ec_velocity)
    if "X_umap" in ec_velocity.obsm:
        scv.pl.velocity_embedding_stream(
            ec_velocity,
            basis="umap",
            color="cell_type_large",
            save="_RNA_velocity.pdf",
        )
    ec_velocity.write_h5ad(OUTPUT_DIR / "endothelial_velocity.h5ad", compression="gzip")




# ==================== 11.主流程 ====================
# 按顺序执行数据处理、注释、EC提取及内皮细胞下游分析。
# 若聚类注释尚未完成，程序会先保存聚类对象并提示填写CLUSTER_TO_CELLTYPE。
def main():
    adata = load_and_cluster()
    adata = annotate_cells(adata)
    _, ec = extract_cell_subsets(adata)
    ec = analyze_endothelial_cells(ec)
    analyze_pathway_and_angiogenesis(ec)
    analyze_tdTomato(ec)
    run_meld(ec)
    analyze_deseq2_results()
    run_velocity(ec)


if __name__ == "__main__":
    main()
