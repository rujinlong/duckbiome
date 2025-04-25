
# DuckBiome project

<!-- badges: start -->

<!-- badges: end -->

This repository contains the complete analysis code used in the
manuscript: “**Bacteria and bacteriophage consortia modulate cecal SCFA
production and host metabolism to enhance feed efficiency in ducks**”.

This study employs a multi-omics approach—integrating 16S amplicon
sequencing, shotgun metagenomics, viromics, liver transcriptomics, and
serum metabolomics—to investigate how gut microbial communities
(bacteria and phages) influence feed efficiency (FE) in ducks. We
identify key microbial taxa, functional pathways, and phage-encoded
auxiliary metabolic genes (AMGs) that correlate with high FE, providing
novel insights into host-microbiome interactions in poultry.

# Dockererized Rstudio server for reproducibility

The dockerfile is in the `docker` folder. A pre-built docker image is
pushed to Docker Hub. Users can pull the image and run the Rstudio
server locally.

``` sh
docker pull jinlongru/duckbiome
docker run -p 80:80 jinlongru/duckbiome
```
