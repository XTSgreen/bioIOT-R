<h1 align="center">bioIOT <span style="font-size:60%">(R)</span></h1>

<p align="center">
  <b>面向单细胞状态转移与拟时序分析的逆最优传输</b>
</p>

<p align="center">
  <a href="https://github.com/XTSgreen/bioIOT-R/actions/workflows/ci.yml"><img src="https://github.com/XTSgreen/bioIOT-R/actions/workflows/ci.yml/badge.svg" alt="R-CMD-check"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="License"></a>
</p>

<p align="center">
  <a href="README.md">English</a> | <b>简体中文</b>
</p>

---

**bioIOT** 是半松弛逆最优传输（IOT）的 R 实现，面向单细胞轨迹分析。给定状态转移特征、源/靶状态质量与观测转移，学习特征权重 θ，使线性代价 `C = -einsum(φ, θ)` 诱导的软边际最优传输计划复现观测数据——并进一步给出状态转移矩阵、随机游走 pseudotime 与 ggplot2 可视化。

求解器在治疗耐药状态转移研究项目中开发并验证，现以通用工具包形式发布。

## 为什么选择 bioIOT？

- **可识别性由构造保证。** 硬边际 OT 下纯列特征不可识别；bioIOT 的 KL 软锚定在保留靶组成约束的同时恢复可识别性。
- **精确隐式梯度。** 前向为纯 base R 的 Anderson 加速不动点求解；反向用隐函数定理——在 unrolled 反传发散之处依然数值稳定。
- **自包含。** 求解器纯 base R 实现；ggplot2 之外无强制运行时依赖。Seurat / SingleCellExperiment 为软门控。
- **端到端。** 细胞嵌入进，转移矩阵 + pseudotime + 图件出。

## 安装

```r
# 从 GitHub 安装
remotes::install_github("XTSgreen/bioIOT-R")
```

<details>
<summary>从本地克隆安装</summary>

```r
install.packages(".", repos = NULL, type = "source")
# 或: R CMD build . && install.packages("bioIOT_0.2.0.tar.gz", repos = NULL)
```

</details>

## 快速上手

```r
library(bioIOT)

# 1) 含真值的可复现合成数据
sim <- simulate_iot_states(K = 6, seed = 1)   # 或 data(demo_iot_states)

# 2) 拟合特征权重：精确隐式梯度 + 两阶段去偏 + 多重启。
#    各场景的状态数 K 可以不同。
fit <- fit_iot(sim$phi, sim$a, sim$b, sim$T_true, n_restart = 2)
fit; summary(fit)

# 3) 轨迹层
Q  <- transition_matrix(fit)                        # (K, K) 转移矩阵
pt <- pseudotime_from_transition(Q, root = "S1")    # 随机游走 pseudotime

# 4) 可视化
plot_transition_heatmap(Q)                          # 带标注热图
plot_transition_flow(Q, sim$embedding)              # CellRank 风格流箭头
plot_theta(fit)                                     # 权重 + 支撑集

# 5) 直接从细胞级对象出发
res <- runIOT(sim$cell_embedding, sim$cell_state,
              from = sim$cell_time == "t0", to = sim$cell_time == "t1",
              root = "S1")
# SingleCellExperiment: runIOT(sce, state_col = "state", time_col = "time",
#                              from = "t0", to = "t1", dimred = "PCA")
# Seurat: runIOT(obj, group.by = "state", split.by = "time",
#                 from = "t0", to = "t1", reduction = "pca")
```

无观测转移 `T_obs` 时，`runIOT()` 以等权特征求解计划；提供 `T_obs`（如谱系/克隆数据）时先拟合权重再求解。

## 实战展示：HSMM 分化时序（真实公开数据）

[`showcase/`](showcase/) 目录在人生骨骼肌成肌细胞分化时序（HSMMSingleCell，271 细胞 × 47k 基因，0/24/48/72 小时；Shin et al. 2015）上完整运行管线，并与 Slingshot 基准对比：

| 拟时序指标（Spearman） | 数值 |
|---|---|
| bioIOT pseudotime ~ 已知时间（细胞级） | 0.251 |
| Slingshot pseudotime ~ 已知时间（细胞级） | 0.253 |
| bioIOT 状态 pseudotime ~ Slingshot 状态 pseudotime | 0.600 |

每状态 30 细胞的采样噪声下，`fit_iot` 恢复真实转移矩阵的精度约为直接使用噪声观测的 3 倍（平均逐行 L1 0.013 vs 0.039；50 次重复）。所有展示数据的来源与许可见
[`showcase/LICENSE_AUDIT.md`](showcase/LICENSE_AUDIT.md)。

## API 总览

| 层 | 函数 |
|---|---|
| 核心求解 | `soft_sinkhorn()`, `row_conditional()`, `make_cost()`, `row_ce_loss()`, `zscore_phi()` |
| 拟合 | `fit_iot()`（含 `print`/`summary` 方法） |
| 轨迹 | `transition_matrix()`, `pseudotime_from_transition()`, `build_state_features()` |
| 单细胞 | `runIOT()` —— matrix / SingleCellExperiment / Seurat |
| 可视化 | `plot_transition_heatmap()`, `plot_transition_flow()`, `plot_theta()`, `plot_pathway_trend()` |
| 演示数据 | `simulate_iot_states()`, `demo_iot_states` |
| 批量队列 | `pathway_markers`, `score_pathways()`, `collapse_probes()`, `gsm_id()` |

## 方法原理

bioIOT 求解

```text
min_P  <C, P> − eps·H(P) + mu·KL(col(P) ‖ b)    s.t.  P·1 = a
```

源端行边际硬、列边际 KL 锚定：

- `mu → ∞` 退化为硬边际 OT（纯列特征不可识别）；
- `mu → 0` 退化为普通行 softmax（失去靶组成锚定）；
- 有限 `mu` 插值两者——论文工作点为 `mu = 0.5, eps = 1.0, lam = 0.05`。

安装后可通过 `browseVignettes("bioIOT")` 查看内置 vignette。

## 测试

```r
library(testthat); library(bioIOT)
test_dir("tests/testthat")   # 仓库根目录
```

## 引用

如果 bioIOT 对你的研究有帮助，请引用：

```bibtex
@misc{dong2026bioiotr,
  author       = {Dong, Han},
  title        = {bioIOT: Inverse Optimal Transport for Single-Cell
                  Trajectory Analysis},
  year         = {2026},
  howpublished = {\url{https://github.com/XTSgreen/bioIOT-R}},
  note         = {R package version 0.2.0}
}
```

## 许可证

[MIT](LICENSE) © 2026 Han Dong (XTSgreen)
