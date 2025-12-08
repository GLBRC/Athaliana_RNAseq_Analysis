# Athaliana_RNAseq_Analysis

Scripts used to analyze RNA-seq data in _A. thaliana_ for Xyloglucan shapes host immunity and pathogen virulence (Bhandari _et al_)

These scripts were used to trim, align, and count the reads on our HTCondor based server cluster:

`bbduk1_submit.txt` = scritp used to trim the reads for quality using BBDUK

`map.sh` and `map.submit` = scripts used to align the RNA-seq reads to the _A. thaliana_ genome using GSNAP

`count.sh` and `count.submit` = scripts used to count the RNA-seq reads that align to _A. thaliana_ genome within defined genes using FeatureCounts

`summarize_counts.sh` and `summarize_counts.submit` = scripts used to summarize mapped and counted reads

`wig.sh` and `wig.submit` = scripts used to generate WIG files for visualization

These scripts were used in R to gernate various graphs and visualizations:

`deseq.R` = R script used to run DESeq2 to identify differentially expressed genes

`xylglucan_plotting_analysis_commands.R` = R script used to analyze the xylglucan data: mean center the data, generate heatmaps, perform and evaulate K-Means clustering, and generate clustering expression plots.

`20250909_XCC_plotting_commands.R` = R script used to analyze the XCC data:  mean center the data, generate heatmaps, perform and evaulate K-Means clustering, and generate clustering expression plots.

`go_enrichment_8categories_logged.R` = R script used to generate GO enrichment from DE gene groups

`plot_of_both_up_down.r` = R script used to generate violin plots

`BP_tsvs_into_heatmap2.r` = R scritp used to generate heatmap of GO BP enrichment

Contact:  dharamch@msu.edu, kmyers2@wisc.edu