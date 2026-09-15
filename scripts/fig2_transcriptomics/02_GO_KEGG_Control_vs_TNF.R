library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

deg_file <- "tables/DEG_Control_vs_TNF_significant.csv"
all_file <- "tables/DEG_Control_vs_TNF_all_results.csv"

deg <- read.csv(deg_file, stringsAsFactors = FALSE)
all_res <- read.csv(all_file, stringsAsFactors = FALSE)

deg_ids <- unique(as.character(deg$GeneID))
universe_ids <- unique(as.character(all_res$GeneID))

# 去掉空值
deg_ids <- deg_ids[!is.na(deg_ids) & deg_ids != ""]
universe_ids <- universe_ids[!is.na(universe_ids) & universe_ids != ""]

# GO 富集分析
ego_bp <- enrichGO(
  gene = deg_ids,
  universe = universe_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_cc <- enrichGO(
  gene = deg_ids,
  universe = universe_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_mf <- enrichGO(
  gene = deg_ids,
  universe = universe_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

write.csv(as.data.frame(ego_bp), "tables/GO_BP_Control_vs_TNF.csv", row.names = FALSE)
write.csv(as.data.frame(ego_cc), "tables/GO_CC_Control_vs_TNF.csv", row.names = FALSE)
write.csv(as.data.frame(ego_mf), "tables/GO_MF_Control_vs_TNF.csv", row.names = FALSE)

# KEGG 富集分析
ekegg <- enrichKEGG(
  gene = deg_ids,
  universe = universe_ids,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

ekegg_readable <- setReadable(
  ekegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

write.csv(as.data.frame(ekegg_readable), "tables/KEGG_Control_vs_TNF.csv", row.names = FALSE)

# 绘图
if (nrow(as.data.frame(ego_bp)) > 0) {
  p1 <- dotplot(ego_bp, showCategory = 20) +
    ggtitle("GO Biological Process: Control vs TNF")
  ggsave("figures/GO_BP_dotplot_Control_vs_TNF.png", p1, width = 9, height = 7, dpi = 300)

  p2 <- barplot(ego_bp, showCategory = 20) +
    ggtitle("GO Biological Process: Control vs TNF")
  ggsave("figures/GO_BP_barplot_Control_vs_TNF.png", p2, width = 9, height = 7, dpi = 300)
}

if (nrow(as.data.frame(ego_cc)) > 0) {
  p3 <- dotplot(ego_cc, showCategory = 20) +
    ggtitle("GO Cellular Component: Control vs TNF")
  ggsave("figures/GO_CC_dotplot_Control_vs_TNF.png", p3, width = 9, height = 7, dpi = 300)
}

if (nrow(as.data.frame(ekegg_readable)) > 0) {
  p4 <- dotplot(ekegg_readable, showCategory = 20) +
    ggtitle("KEGG: Control vs TNF")
  ggsave("figures/KEGG_dotplot_Control_vs_TNF.png", p4, width = 9, height = 7, dpi = 300)

  p5 <- barplot(ekegg_readable, showCategory = 20) +
    ggtitle("KEGG: Control vs TNF")
  ggsave("figures/KEGG_barplot_Control_vs_TNF.png", p5, width = 9, height = 7, dpi = 300)
}

# 提取与屏障、炎症、细胞连接、细胞骨架相关的条目，方便论文撰写
keywords <- c(
  "junction", "adhesion", "barrier", "epithelial",
  "actin", "cytoskeleton", "inflammatory", "inflammation",
  "cytokine", "tumor necrosis factor", "TNF", "NF-kappa", "MAPK"
)

go_bp_df <- as.data.frame(ego_bp)
kegg_df <- as.data.frame(ekegg_readable)

go_bp_focus <- go_bp_df %>%
  filter(str_detect(tolower(Description), paste(tolower(keywords), collapse = "|")))

kegg_focus <- kegg_df %>%
  filter(str_detect(tolower(Description), paste(tolower(keywords), collapse = "|")))

write.csv(go_bp_focus, "tables/GO_BP_focus_barrier_inflammation.csv", row.names = FALSE)
write.csv(kegg_focus, "tables/KEGG_focus_barrier_inflammation.csv", row.names = FALSE)

cat("GO/KEGG enrichment finished.\n")
cat("GO BP terms:", nrow(as.data.frame(ego_bp)), "\n")
cat("GO CC terms:", nrow(as.data.frame(ego_cc)), "\n")
cat("GO MF terms:", nrow(as.data.frame(ego_mf)), "\n")
cat("KEGG pathways:", nrow(as.data.frame(ekegg_readable)), "\n")
cat("Focus GO BP terms:", nrow(go_bp_focus), "\n")
cat("Focus KEGG pathways:", nrow(kegg_focus), "\n")
