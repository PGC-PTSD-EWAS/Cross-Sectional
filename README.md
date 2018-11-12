# PGC EWAS Cross-Sectional

The following scripts are used to run the PGC PTSD EWAS Working Group's cross-sectional analysis based on the methods outlined in [Ratanatharathorn et al. (2017)](https://onlinelibrary.wiley.com/doi/full/10.1002/ajmg.b.32568). There are two steps to the analysis. First, the PGC_PTSD_EWAS_limma.R script is run by each contributing cohort. The results files from this script is sent to the PGC PTSD EWAS Working Group for meta-analysis, where a a series of meta-analysis scripts are performed. Those scripts are:

1. 01_PGC_PTSD_EWAS_loadData.R - this script loads data from each contributing cohort and saves a combined Rdata file with all the data as well as outputting a table with a summary of each cohort's top results. 
