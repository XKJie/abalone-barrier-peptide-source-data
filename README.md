# Submission Figures, Scripts, and Source Data

This folder contains only the scripts and data that were actually used by the article, or are required to reproduce the script-traceable figures in the manuscript.

## Included scope

Only `Figure 1` and `Figure 2` are included here, because these are the figures supported by the audited local scripts.

## Folder layout

- `scripts/`
  Scripts used for the article's script-traceable figures.
- `figures/final_submission/`
  Final submitted `Fig1.jpg` and `Fig2.jpg`.
- `figures/generated_sources/`
  Locally generated panel figures or publication-style figure exports used to assemble the article figures.
- `data/`
  Minimal required source data for reproducing `Figure 1` and `Figure 2`.
- `docs/`
  Panel-level mapping table.

## Figure 1

Used by article:
- final submitted figure: `figures/final_submission/Fig1.jpg`
- generated publication-style figure: `figures/generated_sources/Fig1/Fig3_1_physicochemical_properties.png`
- scripts:
  - `scripts/fig1_peptide_screening/step1_digest.py`
  - `scripts/fig1_peptide_screening/step2_properties.py`
  - `scripts/fig1_peptide_screening/step3_filter.py`
  - `scripts/fig1_peptide_screening/step4_select.py`
  - `scripts/fig1_peptide_screening/step4b_dedup.py`
  - `scripts/fig1_peptide_screening/plot_step4b_publication.py`
- required source data:
  - `data/Fig1/protein.faa`
  - `data/Fig1/step4b_final_dedup.csv`
  - `data/Fig1/panels/A` to `data/Fig1/panels/F`

Files intentionally removed from this package:
- family report and FASTA export, because they are not needed to reproduce the figure shown in the paper.
- extra PDF and SVG figure exports, because the PNG is sufficient for traceability here.

## Figure 2

Used by article:
- final submitted figure: `figures/final_submission/Fig2.jpg`
- generated component figures:
  - `figures/generated_sources/Fig2/A_volcano_Control_vs_TNF.png`
  - `figures/generated_sources/Fig2/B_heatmap_top50_Control_vs_TNF.png`
  - `figures/generated_sources/Fig2/C_barrier_gene_expression_Control_vs_TNF.png`
  - `figures/generated_sources/Fig2/GO_BP_dotplot_publication.png`
  - `figures/generated_sources/Fig2/KEGG_dotplot_publication.png`
- scripts:
  - `scripts/fig2_transcriptomics/01_DESeq2_Control_vs_TNF.R`
  - `scripts/fig2_transcriptomics/02_GO_KEGG_Control_vs_TNF.R`
  - `scripts/fig2_transcriptomics/04_redraw_GO_KEGG_publication.R`
- required source data:
  - `data/Fig2/GSE116936_raw_counts_GRCh38.p13_NCBI.tsv`
  - `data/Fig2/Human.GRCh38.p13.annot.tsv`
  - `data/Fig2/sample_info.csv`
  - `data/Fig2/panels/A` to `data/Fig2/panels/E`

Important note for panel B:
- The exact heatmap matrix is regenerated in script `01_DESeq2_Control_vs_TNF.R` from the raw count matrix plus sample metadata.
- Therefore the retained necessary inputs are the raw count file, annotation file, sample metadata, and the panel B helper files in `data/Fig2/panels/B/`.

## Panel mapping

See:
- `docs/PANEL_DATA_MAPPING.csv`
