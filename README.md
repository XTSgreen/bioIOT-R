# bioIOT (R)

**即插即用**的单细胞时序/拟时序分析包：半松弛逆最优传输（IOT）核心求解器 + 状态转移矩阵 + 随机游走 pseudotime + ggplot2 可视化 + Seurat/SingleCellExperiment 接口 + 批量队列通路工具。核心数值与论文求解器同源。

## 安装

```r
install.packages("packages/bioIOT-R", repos = NULL, type = "source")
# 或 R CMD build packages/bioIOT-R && install.packages("bioIOT_0.2.0.tar.gz", repos = NULL)
```

运行时仅依赖 `ggplot2`；Seurat/SingleCellExperiment 为软门控（装了才启用对应方法）。

## 60 秒上手

```r
library(bioIOT)

# 1) 一行生成可复现合成数据（含真值；或用自带数据 data(demo_iot_states)）
sim <- simulate_iot_states(K = 6, seed = 1)

# 2) 拟合特征权重（隐式微分梯度 + 两阶段去偏 + 多重启）
fit <- fit_iot(sim$phi, sim$a, sim$b, sim$T_true, n_restart = 2, epochs = 150)
fit; summary(fit)

# 3) 状态转移矩阵 + 随机游走 pseudotime
Q <- transition_matrix(fit)
pt <- pseudotime_from_transition(Q, root = "S1")

# 4) 可视化
plot_transition_heatmap(Q)                     # 转移热图
plot_transition_flow(Q, sim$embedding)         # CellRank 风格流箭头
plot_theta(fit)                                # 特征权重条形图

# 5) 直接跑单细胞对象（无 T_obs 用等权解 plan；有 T_obs（如谱系数据）先拟合）
res <- runIOT(sim$cell_embedding, sim$cell_state,
              from = sim$cell_time == "t0", to = sim$cell_time == "t1",
              root = "S1")
res$Q; res$pseudotime
# SingleCellExperiment: runIOT(sce, state_col = "state", time_col = "time", from = "t0", to = "t1")
# Seurat: runIOT(obj, group.by = "state", split.by = "time", from = "t0", to = "t1")
```

底层函数式 API（与论文记号一一对应）：`soft_sinkhorn()`、`make_cost()`、`row_ce_loss()`、`zscore_phi()`、`build_state_features()`。

## 批量队列工具

```r
pw <- score_pathways(expr_gene)                       # 8 维 IOT 通路打分
plot_pathway_trend(pw, time = pseudotime)             # 通路轨迹
expr_gene <- collapse_probes(expr, probe_id, gene_id) # 探针折叠
gsm <- gsm_id(cel_files)                              # CEL 文件名 -> GSM
has_cel_file("GSE62254_RAW.tar")                      # RAW tar 体检
```

## 内容总览

| 层 | 函数 |
|---|---|
| 核心求解 | `soft_sinkhorn`, `row_conditional`, `make_cost`, `row_ce_loss`, `zscore_phi` |
| 拟合 | `fit_iot`（print/summary 方法） |
| 轨迹 | `transition_matrix`, `pseudotime_from_transition`, `build_state_features` |
| 单细胞接口 | `runIOT`（matrix / SingleCellExperiment / Seurat） |
| 可视化 | `plot_transition_heatmap`, `plot_transition_flow`, `plot_theta`, `plot_pathway_trend` |
| 数据 | `simulate_iot_states`, `demo_iot_states` |
| 批量队列 | `pathway_markers`, `score_pathways`, `collapse_probes`, `gsm_id`, `has_cel_file` |

## 质量

- testthat 测试 **84 项全通过**（含隐式梯度 vs 有限差分、θ 恢复、SCE/Seurat 集成）
- `R CMD check --as-cran`：**0 错误 / 1 环境警告（本机缺 qpdf）/ 2 环境性 NOTE**（新提交声明、本机缺 pandoc）
- vignette（Sweave，无需 pandoc）

## 测试与检查

```r
R CMD build packages/bioIOT-R
R CMD check --as-cran bioIOT_0.2.0.tar.gz
# 或直接跑测试
Rscript -e "library(testthat); library(bioIOT); test_dir('packages/bioIOT-R/tests/testthat')"
```
