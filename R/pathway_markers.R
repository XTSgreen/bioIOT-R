# bioIOT pathway marker library
#
# 8 维通路 marker 基因，与论文 IOT 代价特征库同源
# （interpret_drivers.py / 13_bulk_cohort_gex.R）。

#' IOT pathway marker genes (8 pathways)
#'
#' The pathway marker gene library used by the paper's 8-dimensional IOT
#' cost features. Genes are HGNC symbols.
#'
#' @format A named list of 8 character vectors:
#' \describe{
#'   \item{EMT}{VIM, SNAI1, SNAI2, ZEB1, ZEB2, TWIST1, CDH2, FN1, MMP2, MMP9, LOXL2, TGFBI, COL1A1, COL3A1, SPP1}
#'   \item{Angiogenesis}{VEGFA, VEGFB, VEGFC, KDR, FLT1, PGF, ANGPT1, ANGPT2, TEK, PDGFB, PECAM1}
#'   \item{Hypoxia}{HIF1A, LDHA, SLC2A1, CA9, VEGFA, PGK1, ENO1, BNIP3, NDRG1}
#'   \item{Stemness}{LGR5, SOX2, ALDH1A1, CD44, PROM1, KLF4, MYC, NOTCH1, BMI1}
#'   \item{Immune_Cytotoxic}{CD3D, CD3E, CD8A, CXCL9, CXCL10, IFNG, GZMB, PRF1, NKG7}
#'   \item{CellCycle}{MKI67, PCNA, TOP2A, CDK1, CCNB1, CDC20, CCNA2}
#'   \item{TGFb}{TGFB1, TGFBR1, TGFBR2, SMAD2, SMAD3, SMAD4, CTGF, SERPINE1}
#'   \item{Chemokine}{CXCL12, CXCR4, CCL2, CCR2, CCL5, CCL19, CXCL13}
#' }
#' @export
pathway_markers <- list(
  EMT              = c("VIM","SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","CDH2","FN1",
                       "MMP2","MMP9","LOXL2","TGFBI","COL1A1","COL3A1","SPP1"),
  Angiogenesis     = c("VEGFA","VEGFB","VEGFC","KDR","FLT1","PGF","ANGPT1",
                       "ANGPT2","TEK","PDGFB","PECAM1"),
  Hypoxia          = c("HIF1A","LDHA","SLC2A1","CA9","VEGFA","PGK1","ENO1","BNIP3","NDRG1"),
  Stemness         = c("LGR5","SOX2","ALDH1A1","CD44","PROM1","KLF4","MYC","NOTCH1","BMI1"),
  Immune_Cytotoxic = c("CD3D","CD3E","CD8A","CXCL9","CXCL10","IFNG","GZMB","PRF1","NKG7"),
  CellCycle        = c("MKI67","PCNA","TOP2A","CDK1","CCNB1","CDC20","CCNA2"),
  TGFb             = c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","CTGF","SERPINE1"),
  Chemokine        = c("CXCL12","CXCR4","CCL2","CCR2","CCL5","CCL19","CXCL13")
)
