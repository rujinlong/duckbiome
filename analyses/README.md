# Analyses

## Metagenomics

| File | Description |
|------|-------------|
| `01-import.qmd` | Import protein annotation (Bakta, EggNOG) and QC results (CheckM2, GTDB-Tk) into SQLite database |
| `02-humann3.qmd` | Process HUMAnN3 functional profiling output |
| `11-tse.qmd` | TreeSummarizedExperiment (TSE) object creation and data aggregation |
| `12-DAA.qmd` | Differential Abundance Analysis (ANCOMBC) |
| `13-gene_centric.qmd` | Gene-centric analysis with volcano plot visualization |

## RNA-seq

| File | Description |
|------|-------------|
| `21-rnaseq_liver.qmd` | RNA-seq analysis of liver tissue |

## Viromics

| File | Description |
|------|-------------|
| `30-vir_taxonomy.qmd` | Viral taxonomy analysis and annotation import |
| `31-import_virome.qmd` | Virome data import and visualization |

## Integrated analyses

| File | Description |
|------|-------------|
| `60-import.qmd` | Comprehensive data import (DESeq2, phyloseq, Maaslin2) |
| `61-spearman.qmd` | Spearman correlation between FCR and bacterial species |
| `62-FCR.qmd` | Feed Conversion Ratio (FCR) group comparison and statistics |
| `63-16s.qmd` | 16S rRNA alpha diversity analysis (Shannon, Richness) |
| `64-pie.qmd` | EggNOG KEGG Ortholog (KO) pathway visualization |
| `65-network.qmd` | Co-occurrence network analysis (NetCoMi) |
| `66-metabolite.qmd` | Metabolomics data preparation for MetaboAnalyst |
| `67-metabolite_spearman.qmd` | Spearman correlation between metabolites and bacteria |
| `68-average_relative.qmd` | Average relative abundance visualization |

## Misc

| File | Description |
|------|-------------|
| `90-tmp.qmd` | Temporary/scratch notebook |
