---
title: "Drawing coverage plot"
author: "Zaki Fadlullah"
date: "6 April 2016"
output: html_document
---

# transcriptGraph_allGenes

Drawing 5 prime to 3 prime plot for all genes


The workflow invovles three stages
1) Generating appropriate exon level bed file - Using R
  - `GettingExon_region_mouse.R`
2) Calculate individual base coverage - Using bedtools
  - `bedtools_cov.sh`
3) Drawing the actual plot - Using R 
  - `porcess_bedtoolsCov_drawGraph.R`


Final result would look 

![](figures/example_figure.pdf)


