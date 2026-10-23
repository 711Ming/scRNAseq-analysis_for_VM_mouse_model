# 分析流程

```text
10x矩阵
  |
  v
质量控制与标准化
  |
  v
HVG -> PCA -> Harmony -> 邻近图 -> UMAP/Louvain
  |
  v
cluster marker与细胞大类注释
  |
  v
tdTomato状态判定
  |
  v
去除fEC并提取EC
  |
  +--> EC亚群特征、差异表达、细胞周期、动静脉评分
  |
  +--> pseudo-bulk矩阵导出 -> DESeq2结果 -> GO/KEGG/GSEA
  |
  +--> MELD疾病状态评分
  |
  +--> loom匹配 -> RNA velocity
```

细胞类型注释映射需要根据实际marker结果人工确认。
