library(data.table)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)

dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

count_file <- "raw/GSE116936_raw_counts_GRCh38.p13_NCBI.tsv"
annot_file <- "raw/Human.GRCh38.p13.annot.tsv"
sample_file <- "tables/sample_info.csv"

counts_raw <- fread(count_file, data.table = FALSE, check.names = FALSE)
annot <- fread(annot_file, data.table = FALSE, check.names = FALSE)
sample_info <- read.csv(sample_file, stringsAsFactors = FALSE)

colnames(counts_raw)[1] <- "GeneID"

sample_info_sub <- sample_info %>%
  filter(group %in% c("Ctrl", "TNF"))

missing_samples <- setdiff(sample_info_sub$sample, colnames(counts_raw))
if (length(missing_samples) > 0) {
  stop("count matrix missing samples: ", paste(missing_samples, collapse = ", "))
}

count_df <- counts_raw %>%
  select(GeneID, all_of(sample_info_sub$sample)) %>%
  as.data.frame()

rownames(count_df) <- count_df$GeneID
count_df$GeneID <- NULL

count_mat <- round(as.matrix(count_df))
mode(count_mat) <- "integer"

sample_info_sub <- sample_info_sub %>%
  column_to_rownames("sample")

sample_info_sub$group <- factor(sample_info_sub$group, levels = c("Ctrl", "TNF"))

keep <- rowSums(count_mat >= 10) >= 3
count_mat <- count_mat[keep, ]

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = sample_info_sub,
  design = ~ group
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("group", "TNF", "Ctrl"))

res_df <- as.data.frame(res) %>%
  rownames_to_column("GeneID")

annot_sub <- annot %>%
  select(GeneID, Symbol, Description, GeneType)

annot_sub$GeneID <- as.character(annot_sub$GeneID)
res_df$GeneID <- as.character(res_df$GeneID)

res_annot <- res_df %>%
  left_join(annot_sub, by = "GeneID") %>%
  mutate(
    Symbol = ifelse(is.na(Symbol) | Symbol == "", GeneID, Symbol),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 0.5 ~ "Up",
      padj < 0.05 & log2FoldChange < -0.5 ~ "Down",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(padj)

write.csv(
  res_annot,
  "tables/DEG_Control_vs_TNF_all_results.csv",
  row.names = FALSE
)

deg_sig <- res_annot %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5)

write.csv(
  deg_sig,
  "tables/DEG_Control_vs_TNF_significant.csv",
  row.names = FALSE
)

write.table(
  deg_sig$Symbol,
  "tables/DEG_gene_list_for_STRING.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

summary_df <- data.frame(
  comparison = "TNF_vs_Ctrl",
  total_tested_genes = nrow(res_annot),
  significant_DEGs = nrow(deg_sig),
  upregulated = sum(deg_sig$log2FoldChange > 0, na.rm = TRUE),
  downregulated = sum(deg_sig$log2FoldChange < 0, na.rm = TRUE)
)

write.csv(
  summary_df,
  "tables/DEG_summary_Control_vs_TNF.csv",
  row.names = FALSE
)

volcano_df <- res_annot %>%
  mutate(
    neg_log10_padj = -log10(padj),
    label_gene = ifelse(Symbol %in% c("MYLK", "TJP1", "OCLN", "CLDN1", "CLDN2"), Symbol, NA)
  )

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = regulation), alpha = 0.65, size = 1.4) +
  scale_color_manual(
    values = c(
      "Up" = "#D55E00",
      "Down" = "#0072B2",
      "Not significant" = "grey70"
    )
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(aes(label = label_gene), na.rm = TRUE, max.overlaps = 20) +
  theme_bw() +
  labs(
    title = "Control vs TNF",
    x = "log2 Fold Change",
    y = "-log10 adjusted P value"
  )

ggsave(
  "figures/volcano_Control_vs_TNF.png",
  p_volcano,
  width = 7,
  height = 6,
  dpi = 300
)

vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

top50 <- deg_sig %>%
  filter(!is.na(Symbol)) %>%
  slice_head(n = 50)

top50_ids <- top50$GeneID
top50_ids <- top50_ids[top50_ids %in% rownames(vsd_mat)]

if (length(top50_ids) >= 2) {
  mat_top50 <- vsd_mat[top50_ids, ]
  rownames(mat_top50) <- res_annot$Symbol[match(rownames(mat_top50), res_annot$GeneID)]

  ann_col <- data.frame(group = sample_info_sub$group)
  rownames(ann_col) <- rownames(sample_info_sub)

  png("figures/heatmap_top50_Control_vs_TNF.png", width = 1800, height = 2200, res = 300)
  pheatmap(
    mat_top50,
    scale = "row",
    annotation_col = ann_col,
    fontsize_row = 6,
    main = "Top 50 DEGs: Control vs TNF"
  )
  dev.off()
}

barrier_genes <- c("MYLK", "TJP1", "OCLN", "CLDN1", "CLDN2")

barrier_res <- res_annot %>%
  filter(Symbol %in% barrier_genes) %>%
  select(GeneID, Symbol, Description, GeneType, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, regulation)

write.csv(
  barrier_res,
  "tables/barrier_gene_DEG_Control_vs_TNF.csv",
  row.names = FALSE
)

barrier_ids <- barrier_res$GeneID
barrier_ids <- barrier_ids[barrier_ids %in% rownames(vsd_mat)]

if (length(barrier_ids) > 0) {
  barrier_expr <- as.data.frame(t(vsd_mat[barrier_ids, , drop = FALSE]))
  colnames(barrier_expr) <- barrier_res$Symbol[match(barrier_ids, barrier_res$GeneID)]
  barrier_expr$sample <- rownames(barrier_expr)
  barrier_expr$group <- sample_info_sub[barrier_expr$sample, "group"]

  barrier_long <- barrier_expr %>%
    pivot_longer(
      cols = all_of(intersect(barrier_genes, colnames(barrier_expr))),
      names_to = "gene",
      values_to = "vst_expression"
    )

  write.csv(
    barrier_long,
    "tables/barrier_gene_expression_long.csv",
    row.names = FALSE
  )

  p_barrier <- ggplot(barrier_long, aes(x = group, y = vst_expression)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 2) +
    facet_wrap(~ gene, scales = "free_y") +
    theme_bw() +
    labs(
      title = "Barrier-related genes",
      x = "",
      y = "VST expression"
    )

  ggsave(
    "figures/barrier_gene_expression_Control_vs_TNF.png",
    p_barrier,
    width = 8,
    height = 5,
    dpi = 300
  )
}

cat("Finished Control vs TNF differential expression analysis.\n")
cat("Significant DEGs:", nrow(deg_sig), "\n")
cat("Upregulated genes:", summary_df$upregulated, "\n")
cat("Downregulated genes:", summary_df$downregulated, "\n")
cat("Output folders: tables/ and figures/\n")
