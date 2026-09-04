# bioIOT showcase 数据许可审计（LICENSE AUDIT）

审计日期：2026-09-04
审计人：Han Dong (XTSgreen)
范围：bioIOT 展示流程（`packages/showcase/`）所使用的全部外部数据与对照方法。

## 1. 真实数据集：HSMMSingleCell

| 项 | 内容 |
|---|---|
| 数据包 | HSMMSingleCell 1.32.0（Bioconductor Experiment Data Package） |
| 许可证 | **Artistic License 2.0**（包 DESCRIPTION 的 License 字段） |
| 数据内容 | 人骨骼肌成肌细胞分化时间过程 scRNA-seq（HSMM），271 细胞 × ~18k 基因（FPKM），时间点 0 / 24 / 48 / 72 小时 |
| 原始研究 | Shin, J. et al. (2015) *Cell Reports* 11:1432–1444, "Single-Cell RNA-Seq with Waterfall Reveals Molecular Cascades underlying Adult Neurogenesis"（HSMM 数据集随 monocle/HSMMSingleCell 发布） |
| 获取方式 | `BiocManager::install("HSMMSingleCell")`；本仓库不复制、不随包再分发任何原始数据 |
| 再分发条款 | Artistic-2.0 允许使用、修改与再分发，要求保留原始版权与许可声明；本 showcase 仅在用户本地安装该包后进行分析演示，输出图表为派生结果 |
| 结论 | **合规**。GEO（GSE52529）源数据为公开发布数据；通过 Bioconductor 官方渠道获取，许可允许本用途 |

## 2. 基准对照方法：slingshot

| 项 | 内容 |
|---|---|
| 包 | slingshot 2.20.0（Bioconductor） |
| 许可证 | **Artistic License 2.0** |
| 用途 | 仅作为拟时序基准对照（与本包 random-walk pseudotime 相关性比较），不随包再分发 |
| 引用 | Street, K. et al. (2018) *BMC Genomics* 19:477 |
| 结论 | **合规** |

## 3. 合成数据（benchmark_synthetic.R）

由 `bioIOT::simulate_iot_states()` 生成，无外部版权与隐私问题。

## 4. bioIOT 自身

MIT License，Copyright (c) 2026 Han Dong (XTSgreen)。见两包内 LICENSE。

## 5. 遗留事项

- GitHub 仓库 `XTSgreen/bioIOT-R`、`XTSgreen/bioIOT-py` 建立后，`R CMD check` 的 URL 404 NOTE 将消失
- 若未来把任何真实数据子集直接打包（`data-raw/`），必须逐数据集重新审计并在包内保留来源声明
