# Dmel_parallel_temperature_adaptation
This repository contains the raw data as well as the full R script for all the analyses in the following manuscript.

Title: Parallel gene expression evolution in natural and laboratory evolved populations

Authors: Sheng-Kai Hsu1,2,‡, Chaimae Belmouaden1,3,‡, Viola Nolte1 and Christian Schlötterer1,*

1 Institut für Populationsgenetik, Vetmeduni Vienna, Vienna, Austria.
2 Vienna Graduate School of Population Genetics, Vetmeduni Vienna, Vienna, Austria.
3 Current affiliation: Faculty of fundamental and applied sciences of Poitiers, France.
‡ The authors contribute equally to the manuscript.
* Correspondence: christian.schloetterer@vetmeduni.ac.at; Tel.: +43-1-25077-4300.

Abstract

Ecological adaptation is frequently inferred by the comparison of natural populations from different environments. Nevertheless, the inference of the selective forces suffers the challenge that many environmental factors covary. With well-controlled environmental conditions, experimental evolution provides a powerful approach to complement the analysis of natural populations. On the other hand, it is apparent that laboratory conditions differ in many ways from natural environments, which raises the question to what extent selection responses in experimental evolution studies can inform us about adaptation processes in the wild. In this study, we compared the expression profiles of replicated Drosophila melanogaster populations which have been exposed to two distinct temperature regimes (18/28 ºC and 10/20 ºC) in the laboratory for more than 80 generations. Using gene-wise differential expression analysis and co-expression network analysis, we identified 541 genes and three co-regulated gene modules that evolved in the same direction in both temperature regimes, and most of these changes probably reflect an adaptation to the space constrain or diurnal temperature fluctuation that is common in both selection regimes. 203 genes and seven modules evolved temperature-specific expression changes. Remarkably, we detected a significant overlap of these temperature-adaptive genes/modules from experimental evolution with temperature-adaptive genes inferred from natural Drosophila populations covering two different temperature clines. We conclude that well-designed experimental evolution studies are a powerful tool to dissect evolutionary responses. 

# Brief Readme for the files
The Rscript is named "Rscript_analysis_published.R"
"count_table_published.csv" contains the read counts on each gene for each sample.
"mel-cline-gene-count.deseq2.SH-NH.exp.chr" is the DE analysis between Northen and Southern North American populations from Zhao et al., 2015.
"hutter_2008_Dmel_22.csv" is the table containing the expression difference between African and Europian populations from Hutter et al., 2008.
"Afr_EU.txt" are the genes showing higher expression in African populations and "EU_Afr.txt" are the genes showing higher expression in European populations.
