########################################################################################
# PGC Methylation Meta-Analysis
########################################################################################

rm(list=ls())

setwd("") # Directory to save results files
homeDir<-"" # Directory with the results from each site

########################################################################################
# Step 1: Load Summary Statistics 
########################################################################################

# Step 1A: Load DNHS Results
load("DNHS100_limma_results_08.09.16.Rdata")
DNHS.coef<-fit.coef
rm(fit.coef)
results<-results[order(results[, "P.Value"]), ]
rank<-c(1:nrow(results))
results<-cbind(results, rank)
DNHS.results<-results
rm(results, rank)
DNHS.ebayes<-fit2.ebayes
rm(fit2.ebayes)

# Step 1B: Load GTP Results
load("GTP_limma_results_02.26.17.Rdata")
GTP.coef<-fit.coef
rm(fit.coef)
results<-results[order(results[, "P.Value"]), ]
GTP.results<-results
rm(results)
rank<-c(1:nrow(GTP.results))
GTP.results<-cbind(GTP.results, rank)
rm(rank)
GTP.ebayes<-fit2.ebayes
rm(fit2.ebayes)

# Step 1C: Load MRS Results
load("MRS_limma_results_oct12_2015.Rdata")
MRS.coef<-fit.coef
rm(fit.coef)
results<-results[order(results[, "P.Value"]), ]
rank<-c(1:nrow(results))
results<-cbind(results, rank)
MRS.results<-results
rm(results, rank)
MRS.ebayes<-fit2.ebayes
rm(fit2.ebayes)

# Step 1D: Load Boston-VA Results
load("BostonVA_Dec_21_2015_limma_results.Rdata")
VA.coef<-fit.coef
rm(fit.coef)
results<-results[order(results[, "P.Value"]), ]
rank<-c(1:nrow(results))
results<-cbind(results, rank)
VA.results<-results
rm(results, rank)
VA.ebayes<-fit2.ebayes
rm(fit2.ebayes)

# Step 1E: Load PRISMO Data
load("utrecht3_limma_results.Rdata")
PRISMO.coef<-fit.coef
rm(fit.coef)
results<-results[order(results[, "P.Value"]), ]
rank<-c(1:nrow(results))
results<-cbind(results, rank)
PRISMO.results<-results
rm(results, rank)
PRISMO.ebayes<-fit2.ebayes
rm(fit2.ebayes)

# Step 1F: Load WTC Data
load("WTC180_meta_limma_results_08.09.16.Rdata")
WTC.coef<-fit.coef
rm(fit.coef)
results<-results[order(results[, "P.Value"]), ]
rank<-c(1:nrow(results))
results<-cbind(results, rank)
WTC.results<-results
rm(results, rank)
WTC.ebayes<-fit2.ebayes
rm(fit2.ebayes)

# Step 1G: Load Duke VA Data
load("Duke_limma_results_103015.Rdata")
DUKE.coef<-fit.coef
rm(fit.coef)
results<-results[order(results[, "P.Value"]), ]
rank<-c(1:nrow(results))
results<-cbind(results, rank)
DUKE.results<-results
rm(results, rank)
DUKE.ebayes<-fit2.ebayes
rm(fit2.ebayes)

# Step 1H: Load Army STARRS Data
load("STARRS_limma_results_mar7_2016.Rdata")
AS.coef<-fit.coef
rm(fit.coef)
results<-results[order(results[, "P.Value"]), ]
rank<-c(1:nrow(results))
results<-cbind(results, rank)
AS.results<-results
rm(results, rank)
AS.ebayes<-fit2.ebayes
rm(fit2.ebayes)

# Step 1I: Load MIRECC-AA Data
load("Duke_AAset1_limma_results_031616.Rdata")
MIR.coef<-fit.coef
rm(fit.coef)
results<-results[order(results[, "P.Value"]), ]
rank<-c(1:nrow(results))
results<-cbind(results, rank)
MIR.results<-results
rm(results, rank)
MIR.ebayes<-fit2.ebayes
rm(fit2.ebayes)

setwd(homeDir)
save.image("PGC_EWAS_combinedData.Rdata")
rm(list=ls())

########################################################################################
# # Step 2: QQ Plots of Individual Data Sources
########################################################################################

load("PGC_EWAS_combinedData.Rdata")

studies<-c("AS", "DNHS", "DUKE", "GTP", "MRS", "MIR", "PRISMO", "VA", "WTC")
tab<-matrix(nrow=length(studies), ncol=9)
colnames(tab)<-c("Study", "Sites", "Min N", "Max N", "Lambda", "FDR", "p<5x10^-5", 
  "p<5x10^-6", "p<5x10^-7")

for(ii in 1:length(studies)){
  tab[ii, "Study"]<-studies[ii]
  res<-get(paste(studies[ii], ".results", sep=""))
  res<-res[!is.na(res[, "logFC"]),]
  tab[ii, "p<5x10^-5"]<-sum(res[, "P.Value"]<(5*10^(-5))) 
  tab[ii, "p<5x10^-6"]<-sum(res[, "P.Value"]<(5*10^(-6))) 
  tab[ii, "p<5x10^-7"]<-sum(res[, "P.Value"]<(5*10^(-7))) 
  tab[ii, "FDR"]<-sum(res[, "adj.P.Val"]<=0.05)
  tab[ii, "Sites"]<-nrow(res) 
  coef<-get(paste(studies[ii], ".coef", sep=""))
  tab[ii, "Min N"]<-range(coef[, "N.subjects"])[1] 
  tab[ii, "Max N"]<-range(coef[, "N.subjects"])[2] 
  p<-res[, "P.Value"]
  chisq<-qchisq(1-p,1)
  tab[ii, "Lambda"]<-round(median(chisq)/qchisq(0.5,1),3)
  rm(p, chisq, coef, res)
}

write.csv(tab, "PGC_EWAS_individualCohort_results.csv")
