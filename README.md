# 小鼠SWS单细胞转录组分析

本仓库整理小鼠Sturge–Weber综合征（SWS）相关单细胞转录组分析流程。代码按原始数据处理、质量控制、聚类与细胞注释、内皮细胞（EC）提取及亚群分析、tdTomato分析、功能富集和RNA velocity的顺序组织。

> 本仓库不包含原始测序数据、`.h5ad`文件、loom文件或差异分析结果。运行前请将数据放入相应目录，并检查脚本中的样本及参数设置。

## 目录结构

```text
mouse-sws-scrnaseq/
├── README.md
├── requirements.txt
├── .gitignore
├── scripts/
│   └── mouse_sws_scrnaseq.py
├── data/
│   ├── matrix/
│   │   └── README.md
│   └── loom/
│       └── README.md
└── results/
    └── .gitkeep
```

## 分析内容

1. 读取10x Genomics表达矩阵并整合样本。
2. 质量控制、标准化、高变基因筛选、PCA、Harmony批次校正、UMAP及Louvain聚类。
3. 导出cluster marker并结合标志基因完成细胞大类注释。
4. 识别tdTomato表达状态并统计不同样本和细胞群中的分布。
5. 去除fEC并提取aEC、vEC、cEC和dEC等内皮细胞亚群。
6. 进行内皮细胞标志基因、dEC差异表达、细胞周期及动静脉特征评分分析。
7. 导出pseudo-bulk计数矩阵及样本信息，供DESeq2分析使用。
8. 对差异表达基因进行GO/KEGG富集及GSEA分析。
9. 进行MELD疾病状态评分。
10. 基于loom文件开展RNA velocity分析。

## 环境配置

建议使用Python 3.10的独立环境。

```bash
python -m venv .venv
```

激活环境后安装依赖：

```bash
pip install -r requirements.txt
```

部分依赖（尤其是`meld`及其底层依赖）可能因操作系统和Python版本而需要单独处理。若不运行MELD或RNA velocity相关分析，可根据实际需要调整相应依赖。

## 输入数据

### 1.10x表达矩阵

每个样本对应一个目录，目录内应包含10x Genomics格式的矩阵文件，例如`matrix.mtx`、`features.tsv`和`barcodes.tsv`，也可以使用相应的压缩文件。

默认样本目录如下：

```text
data/matrix/
├── CM12_1/
├── CM12_2/
├── CM12_3/
├── CM12_4/
├── CM28_1/
├── CM28_2/
└── WT28/
```

样本名称、分组和目录映射在` scripts/mouse_sws_scrnaseq.py`开头的`SAMPLES`变量中设置。请根据实际数据修改。

### 2.loom文件

RNA velocity分析使用的loom文件放在：

```text
data/loom/
```

文件名应与脚本中的`LOOM_SAMPLES`列表一致，例如`CM12_1.loom`。请核对样本条码与10x表达矩阵中的细胞条码能否匹配。

### 3.细胞类型注释

运行全数据集聚类后，脚本会导出`results/cluster_marker_genes.csv`及候选marker可视化。请根据marker表达结果确认每个Louvain cluster对应的细胞类型，并在`CLUSTER_TO_CELLTYPE`中填写映射。例如：

```python
CLUSTER_TO_CELLTYPE = {
    "0": "EC",
    "1": "PC",
    "2": "SMC",
}
```

以上仅为格式示例，不代表实际数据中的注释结果。内皮细胞亚群标签应与`EC_TYPES`保持一致。

## 运行方法

在仓库根目录下运行：

```bash
python scripts/mouse_sws_scrnaseq.py
```

首次运行会先执行原始数据处理、质量控制、降维聚类并导出marker。若尚未填写`CLUSTER_TO_CELLTYPE`，程序会提示补充细胞注释映射；完成后再次运行，以继续完成细胞注释、EC提取和下游分析。

## 主要输出

分析结果写入`results/`，主要包括：

- 质量控制统计表及QC图；
- 聚类结果、UMAP图及cluster marker表；
- 细胞类型和EC亚群的样本计数表；
- tdTomato阳性细胞统计及差异表达结果；
- EC亚群差异基因和功能富集结果；
- pseudo-bulk计数矩阵及样本信息；
- MELD评分结果；
- RNA velocity对象及流线图。

图形文件存放在`results/figures/`。

## 复现与解释注意事项

- 代码中的样本列表及质量控制参数应与最终实际分析记录保持一致。
- 当前默认线粒体基因比例阈值为15%；若最终分析使用其他阈值，应同步修改代码和方法学描述。
- `CLUSTER_TO_CELLTYPE`必须根据实际marker结果人工确认，不能直接使用示例映射。
- pseudo-bulk分析应使用原始整数计数。若输入对象缺少原始计数层，代码不会用log-normalized表达量替代原始计数。
- 统计推断应以独立生物学样本为单位考虑重复结构；不能仅依据单个细胞数量判断生物学重复数。
- 本代码基于已提供的分析Notebook整理。语法检查不等同于使用原始数据完成了端到端运行验证。

## 数据与隐私

原始测序数据、临床信息、个人路径、较大的中间文件和分析结果默认不应直接提交到公开仓库。请确认数据共享权限后再上传，并通过`.gitignore`排除不应公开的文件。
