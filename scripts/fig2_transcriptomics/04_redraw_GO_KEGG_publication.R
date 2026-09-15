#!/usr/bin/env Rscript

library(tidyverse)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

parse_gene_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(z) {
    as.numeric(z[1]) / as.numeric(z[2])
  })
}

wrap_label <- function(x, width = 36) {
  stringr::str_wrap(x, width = width)
}

make_enrichment_dotplot <- function(
    infile,
    outfile_prefix,
    plot_title,
    top_n = 12,
    label_width = 36
) {
  df <- read.csv(infile, stringsAsFactors = FALSE)

  if (nrow(df) == 0) {
    message("No enrichment terms found in: ", infile)
    return(NULL)
  }

  df_plot <- df %>%
    filter(!is.na(p.adjust)) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      GeneRatioNumeric = parse_gene_ratio(GeneRatio),
      DescriptionWrapped = wrap_label(Description, width = label_width),
      DescriptionWrapped = factor(
        DescriptionWrapped,
        levels = rev(DescriptionWrapped)
      ),
      negLog10Padj = -log10(p.adjust)
    )

  p <- ggplot(
    df_plot,
    aes(
      x = GeneRatioNumeric,
      y = DescriptionWrapped,
      size = Count,
      color = negLog10Padj
    )
  ) +
    geom_point(alpha = 0.90) +
    scale_color_gradient(
      low = "#5B8DB8",
      high = "#C44E52",
      name = expression(-log[10]("adjusted P"))
    ) +
    scale_size_continuous(
      range = c(3.0, 7.5),
      name = "Gene count"
    ) +
    scale_x_continuous(
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      title = plot_title,
      x = "Gene ratio",
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 13
      ),
      axis.text.y = element_text(
        size = 10,
        color = "black",
        lineheight = 0.9
      ),
      axis.text.x = element_text(
        size = 10,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 11,
        color = "black",
        margin = margin(t = 8)
      ),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.key.height = unit(0.45, "cm"),
      legend.key.width = unit(0.35, "cm"),
      plot.margin = margin(8, 12, 8, 8)
    )

  ggsave(
    paste0(outfile_prefix, ".png"),
    p,
    width = 7.2,
    height = 5.4,
    dpi = 600
  )

  ggsave(
    paste0(outfile_prefix, ".pdf"),
    p,
    width = 7.2,
    height = 5.4
  )

  message("Saved: ", outfile_prefix, ".png / .pdf /")
}

make_enrichment_dotplot(
  infile = "tables/GO_BP_Control_vs_TNF.csv",
  outfile_prefix = "figures/GO_BP_dotplot_publication",
  plot_title = "GO biological process",
  top_n = 12,
  label_width = 36
)

make_enrichment_dotplot(
  infile = "tables/KEGG_Control_vs_TNF.csv",
  outfile_prefix = "figures/KEGG_dotplot_publication",
  plot_title = "KEGG pathways",
  top_n = 12,
  label_width = 36
)

message("Publication-style GO and KEGG dotplots finished.")
