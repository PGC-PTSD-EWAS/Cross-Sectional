########################################################################################
# PGC EWAS Meta-Analysis: Summarizing Results
########################################################################################

rm(list=ls())
library(rmeta)
library(ChAMP)
library(forestplot)
library(qqman)
library(ggplot2)

setwd("") # Directory with the results from scripts 1-2

########################################################################################
# Step 1: Loading Data
########################################################################################

# Combining results files
load("PGC_EWAS_nonSmoke_inverseNorm_intersectSites.Rdata")
z.intersect<-z.combined
rm(z.combined)
load("PGC_EWAS_nonSmoke_inverseNorm_remainingSites.Rdata")
z.remain<-z.combined
rm(z.combined)
z.stats<-c(z.intersect, z.remain)
rm(z.intersect, z.remain)

all<-pval.combined.marot
rm(pval.combined.marot)
colnames(all)<-c("CpG", "p")

# Replace this with the cohorts included in the analysis
Studies<-rep(c("DNHS, GTP, MRS, VA, PRISMO, WTC, DUKE, AS, MIR, TRUST"), nrow(all)) 
all<-cbind(all, Studies)
all<-all[, colnames(results)]
all(colnames(all)==colnames(results))
results<-rbind(all, results)

# Removing CpG sites not in at least two studies
results<-results[results[, "Studies"]!="DNHS",]
results<-results[results[, "Studies"]!="GTP",]
results<-results[results[, "Studies"]!="MRS",]
results<-results[results[, "Studies"]!="VA",]
results<-results[results[, "Studies"]!="PRISMO",]
results<-results[results[, "Studies"]!="WTC",]
results<-results[results[, "Studies"]!="DUKE",]
results<-results[results[, "Studies"]!="AS",]
results<-results[results[, "Studies"]!="MIR",]
results<-results[results[, "Studies"]!="TRUST",]

results<-results[order(as.numeric(results[, "p"])), ]
FDR<-p.adjust(as.numeric(results[, "p"]), method="fdr", n=nrow(results))
results<-cbind(results, FDR)

pvalue <- as.numeric(results[, "p"])
chisq <- qchisq(1-pvalue,1)
lambda = median(chisq)/qchisq(0.5,1) 

z<-data.frame(z.stats)
z$norm<-rnorm(n=nrow(z), mean=0, sd=1)

png("PGC_EWAS_inverseNorm_histogram_wINTRuST.png", width=600, height=600, units="px", bg="white")
ggplot(z, aes(x=z.stats))+geom_histogram(aes(y = ..density..))+geom_density(data=z, aes(x=norm))+
  ggtitle(paste("PGC EWAS w/ INTRuST, l = ", round(lambda, 3), sep=""))+xlab("test-statistics")+ylab("Density")+
  theme_bw()+theme(text=element_text(size=20))
dev.off()

png("PGC_EWAS_inverseNorm_histogram_wINTRuST_publication.png", width=1200, height=600, units="px", bg="white")
ggplot(z, aes(x=z.stats))+geom_histogram(aes(y = ..density..))+geom_density(data=z, aes(x=norm))+
  ggtitle("")+xlab("test-statistics")+ylab("Density")+
  theme_bw()+theme(text=element_text(size=24))
dev.off()


save(results, file="PGC_EWAS_inverseNorm_allResults.Rdata")

# Summary info
tab<-matrix(nrow=9, ncol=2)
colnames(tab)<-c("Parameter", "Value")
tab[, "Parameter"]<-c("CpG sites", "N sites in 3 studies", "N sites in 2 studies",
                      "N sites in 1 study", "N sites p<5x10^-5", "N sites p<5x10^-6",
                      "N sites p<5x10^-7", "N sites FDR < 0.05", "lambda")
rownames(tab)<-tab[,"Parameter"]
pvalue <- as.numeric(results[, "p"])
chisq <- qchisq(1-pvalue,1)
lambda = median(chisq)/qchisq(0.5,1)
tab["lambda", "Value"]<-round(lambda,6)

tab["N sites FDR < 0.05", "Value"]<-sum(FDR<=0.05)
tab["N sites p<5x10^-5", "Value"]<-sum(pvalue<=5*10^-5)
tab["N sites p<5x10^-6", "Value"]<-sum(pvalue<=5*10^-6) # 7
tab["N sites p<5x10^-7", "Value"]<-sum(pvalue<=5*10^-7) # 3

nStudies<-unlist(lapply(strsplit(results[, "Studies"], ", "), function(x) length(x)))
tab["N sites in 3 studies", "Value"]<-sum(nStudies==3)
tab["N sites in 2 studies", "Value"]<-sum(nStudies==2)
tab["N sites in 1 study", "Value"]<-sum(nStudies==1)
tab["CpG sites", "Value"]<-nrow(results)
rownames(tab)<-c(1:nrow(tab))

write.csv(tab, file="PGC_EWAS_inverseNorm_allResults_summary.csv")
rm(list=ls()[-match("results", ls())])

########################################################################################
# Step 2: QQ plot
########################################################################################

ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=2, cex.lab=4.5, cex.axis=4, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}
pvalue <- as.numeric(results[, "p"])

png("PGC_EWAS_inverseNorm_QQplot_wINTRuST.png", width=1600, height=800, units="px", bg="white")
par(mar=c(10,10,4,2), mgp=c(6,2,0))
ggd.qqplot(pvalue, main = paste(""))
dev.off()

rm(list=ls()[-match("results", ls())])

########################################################################################
# Step 3: Manhattan plot
########################################################################################

rm(list=ls())

load("PGC_EWAS_inverseNorm_allResults.Rdata")
rownames(results)<-results[, "CpG"]
data(probe.features)
results<-results[which(rownames(results)%in%rownames(probe.features)),]
results<-results[-grep("rs", rownames(results)),]
probe.features$CHR<-as.character(probe.features$CHR)
probe.features<-probe.features[rownames(results),]
probe.features$CHR[probe.features$CHR=="X"]<-23
probe.features$CHR[probe.features$CHR=="Y"]<-24
probe.features$CHR<-as.numeric(probe.features$CHR)
probe.features<-as.matrix(probe.features[, c("CHR", "MAPINFO", "gene")])
all(rownames(results)==rownames(probe.features))
results<-cbind(results, probe.features)

head(results)
df<-data.frame(P=as.numeric(results[, "p"]), 
               FDR=as.numeric(results[, "FDR"]),
               CHR=as.numeric(results[, "CHR"]), 
               BP=as.numeric(results[, "MAPINFO"]),
               SNP=results[, "gene"], stringsAsFactors=F)


df["cg17284326", "SNP"]<-"MIR3170"
df["cg26801037", "SNP"]<-"AC011899"
df["cg18217048", "SNP"]<-"LINC00599"
cutpoint<--log10(1.210462*10^-06)
cut2<--log10(2.406347*10^-05)
labs<-c(1:22, "X", "Y")

opar<-par()

t<-results[results[, "FDR"]<=0.05,]
t<-probe.features[rownames(t),]
write.csv(t, "PGC_EWAS_topresults.csv")

png("PGC_EWAS_manhattan_wINTRuST.png",height=600, width=1200, units="px")
par(mar=c(5.1,6.1, 4.1, 2.1))
manhattan(df, suggestiveline=FALSE, genomewideline=cutpoint, ylim=c(0,12),
          chrlabs=labs, cex=2, cex.axis=2, cex.lab=3, col=c("blue4", "red4"))
dev.off()

postscript("PGC_EWAS_manhattan_wINTRuST.eps", horizontal=FALSE, width = 6, height = 3, 
           onefile = FALSE, paper = "special", colormodel = "cmyk")
#par(mar=c(5.1,6.1, 4.1, 2.1))
par(opar)
manhattan(df, suggestiveline=FALSE, genomewideline=FALSE, ylim=c(0,12),
          chrlabs=labs, cex=0.5, cex.axis=0.5, cex.lab=1)#, col=c("blue4", "red4"))
abline(h=cutpoint, lty=2)
abline(h=cut2, lty=3)
dev.off()

png("PGC_EWAS_manhattan_wINTRuST_BW.png", width = 600, height = 300, units="px")
manhattan(df, suggestiveline=FALSE, genomewideline=FALSE, ylim=c(0,12),
          chrlabs=labs, cex=0.5, cex.axis=0.75, cex.lab=1)#, col=c("blue4", "red4"))

abline(h=cutpoint, lty=2)
abline(h=cut2, lty=3)
dev.off()

png("PGC_EWAS_manhattan_wINTRuST_BW_ANNOTATED.png", width = 600, height = 300, units="px")
manhattan(df, suggestiveline=FALSE, genomewideline=FALSE, ylim=c(0,12),
          annotatePval = 1.210462*10^-06, annotateTop = F,
          chrlabs=labs)#, cex=0.5, cex.axis=0.75, cex.lab=1)#, col=c("blue4", "red4"))

abline(h=cutpoint, lty=2)
abline(h=cut2, lty=3)
dev.off()

rm(list=ls())

########################################################################################
# Step 4: Top Results Table
########################################################################################

load("PGC_EWAS_inverseNorm_allResults.Rdata")
rownames(results)<-results[, "CpG"]
cpgs<-rownames(results[as.numeric(results[, "p"])<=5*10^-5, ]) 
topResults<-c("cg05575921", "cg21161138", "cg25648203", "cg26703534", 
  "cg25415650", "cg17284326")
topResults<-topResults[!topResults%in%cpgs]

# From the CHARGE study, load the most siginificant smoking and former smoking sites
sigSmoking<-read.csv("Most_significant_smoking_cpg.csv",
                     stringsAsFactors=F)
sigSmoking<-sigSmoking$Probe.ID
sigSmoking<-sigSmoking[sigSmoking%in%rownames(results)]
sigSmoking<-sigSmoking[!sigSmoking%in%cpgs]

formerSmoking<-read.csv("former_smoking_sites.csv",
                        stringsAsFactors=F)
formerSmoking<-formerSmoking$Probe.ID
formerSmoking<-formerSmoking[formerSmoking%in%rownames(results)]
formerSmoking<-formerSmoking[!formerSmoking%in%cpgs]
formerSmoking<-formerSmoking[!formerSmoking%in%sigSmoking]

results<-results[c(cpgs, topResults, sigSmoking, formerSmoking), ] 
results<-data.frame(results, stringsAsFactors=F)
results$p<-as.numeric(results$p)
results$FDR<-as.numeric(results$FDR)
sites<-rownames(results)

load("PGC_EWAS_DataPrep_nonSmoke.Rdata")
rm(list=ls()[grep("ebayes", ls())])

colnames(DNHS.coef)
DNHS.coef<-DNHS.coef[rownames(DNHS.coef)%in%sites,c("PTSDpm", "N.subjects")]
DNHS.results<-DNHS.results[rownames(DNHS.coef),]
colnames(DNHS.coef)<-c("PTSD", "N.subjects") # need consistent column headings
all(rownames(DNHS.coef)==rownames(DNHS.results))
DNHS<-data.frame(DNHS.coef, DNHS.results)
str(DNHS) # should all be numbers not factors
DNHS$s.e.<-DNHS$PTSD/DNHS$t # calculate the standard error
DNHS$weight<-1/(DNHS$s.e.^2) # calculate weight
rm(DNHS.results, DNHS.coef, DNHS.oneSided)

colnames(GTP.coef)
GTP.coef<-GTP.coef[rownames(GTP.coef)%in%sites,c("PTSDcurr", "N.subjects")]
GTP.results<-GTP.results[rownames(GTP.coef),]
colnames(GTP.coef)<-c("PTSD", "N.subjects")
all(rownames(GTP.coef)==rownames(GTP.results))
GTP<-data.frame(GTP.coef, GTP.results)
str(GTP)
GTP$s.e.<-GTP$PTSD/GTP$t # calculate the standard error
GTP$weight<-1/(GTP$s.e.^2) # calculate weight
rm(GTP.results, GTP.coef, GTP.oneSided)


colnames(MRS.coef)
MRS.coef<-MRS.coef[rownames(MRS.coef)%in%sites,c("PTSDbroad", "N.subjects")]
MRS.results<-MRS.results[rownames(MRS.coef),]
colnames(MRS.coef)<-c("PTSD", "N.subjects")
all(rownames(MRS.coef)==rownames(MRS.results))
MRS<-data.frame(MRS.coef, MRS.results)
str(MRS)
MRS$s.e.<-MRS$PTSD/MRS$t # calculate the standard error
MRS$weight<-1/(MRS$s.e.^2) # calculate weight
rm(MRS.results, MRS.coef, MRS.oneSided)

colnames(VA.coef)
VA.coef<-VA.coef[rownames(VA.coef)%in%sites,c("PTSD_C", "N.subjects")]
VA.results<-VA.results[rownames(VA.coef),]
colnames(VA.coef)<-c("PTSD", "N.subjects")
all(rownames(VA.coef)==rownames(VA.results))
VA<-data.frame(VA.coef, VA.results)
str(VA)
VA$s.e.<-VA$PTSD/VA$t # calculate the standard error
VA$weight<-1/(VA$s.e.^2) # calculate weight
rm(VA.results, VA.coef, VA.oneSided)

colnames(PRISMO.coef)
PRISMO.coef<-PRISMO.coef[rownames(PRISMO.coef)%in%sites,c("ptssTRUE", "N.subjects")]
PRISMO.results<-PRISMO.results[rownames(PRISMO.coef),]
colnames(PRISMO.coef)<-c("PTSD", "N.subjects")
all(rownames(PRISMO.coef)==rownames(PRISMO.results))
PRISMO<-data.frame(PRISMO.coef, PRISMO.results)
str(PRISMO)
PRISMO$s.e.<-PRISMO$PTSD/PRISMO$t # calculate the standard error
PRISMO$weight<-1/(PRISMO$s.e.^2) # calculate weight
rm(PRISMO.results, PRISMO.coef, PRISMO.oneSided)

colnames(WTC.coef)
WTC.coef<-WTC.coef[rownames(WTC.coef)%in%sites,c("PTSD", "N.subjects")]
WTC.results<-WTC.results[rownames(WTC.coef),]
colnames(WTC.coef)<-c("PTSD", "N.subjects")
all(rownames(WTC.coef)==rownames(WTC.results))
WTC<-data.frame(WTC.coef, WTC.results)
str(WTC)
WTC$s.e.<-WTC$PTSD/WTC$t # calculate the standard error
WTC$weight<-1/(WTC$s.e.^2) # calculate weight
rm(WTC.results, WTC.coef, WTC.oneSided)

colnames(DUKE.coef)
DUKE.coef<-DUKE.coef[rownames(DUKE.coef)%in%sites,c("currptsd", "N.subjects")]
DUKE.results<-DUKE.results[rownames(DUKE.coef),]
colnames(DUKE.coef)<-c("PTSD", "N.subjects")
all(rownames(DUKE.coef)==rownames(DUKE.results))
DUKE<-data.frame(DUKE.coef, DUKE.results)
str(DUKE)
DUKE$s.e.<-DUKE$PTSD/DUKE$t # calculate the standard error
DUKE$weight<-1/(DUKE$s.e.^2) # calculate weight
rm(DUKE.results, DUKE.coef, DUKE.oneSided)

colnames(AS.coef)
AS.coef<-AS.coef[rownames(AS.coef)%in%sites,c("d_pts30_t2", "N.subjects")]
AS.results<-AS.results[rownames(AS.coef),]
colnames(AS.coef)<-c("PTSD", "N.subjects")
all(rownames(AS.coef)==rownames(AS.results))
AS<-data.frame(AS.coef, AS.results)
str(AS)
AS$s.e.<-AS$PTSD/AS$t # calculate the standard error
AS$weight<-1/(AS$s.e.^2) # calculate weight
rm(AS.results, AS.coef, AS.oneSided)

colnames(MIR.coef)
MIR.coef<-MIR.coef[rownames(MIR.coef)%in%sites,c("curr_ptsd", "N.subjects")]
MIR.results<-MIR.results[rownames(MIR.coef),]
colnames(MIR.coef)<-c("PTSD", "N.subjects")
all(rownames(MIR.coef)==rownames(MIR.results))
MIR<-data.frame(MIR.coef, MIR.results)
str(MIR)
MIR$s.e.<-MIR$PTSD/MIR$t # calculate the standard error
MIR$weight<-1/(MIR$s.e.^2) # calculate weight
rm(MIR.results, MIR.coef, MIR.oneSided)

colnames(TRUST.coef)
TRUST.coef<-TRUST.coef[rownames(TRUST.coef)%in%sites,c("PTSD", "N.subjects")]
TRUST.results<-TRUST.results[rownames(TRUST.coef),]
colnames(TRUST.coef)<-c("PTSD", "N.subjects")
all(rownames(TRUST.coef)==rownames(TRUST.results))
TRUST<-data.frame(TRUST.coef, TRUST.results)
str(TRUST)
TRUST$s.e.<-TRUST$PTSD/TRUST$t # calculate the standard error
TRUST$weight<-1/(TRUST$s.e.^2) # calculate weight
rm(TRUST.results, TRUST.coef, TRUST.oneSided)

head(results[,"Studies"])
studies<-unlist(strsplit(results[1, "Studies"], ","))
studies<-gsub(" ", "", studies)

for(ii in 1:length(studies)){
  mat<-matrix(nrow=nrow(results), ncol=3)
  rownames(mat)<-rownames(results)
  colnames(mat)<-paste(studies[ii], c(".beta", ".se", ".N"), sep="")
  results<-cbind(results, mat)
  rm(mat)
}

results$variance<-results$beta<-NA

for(ii in 1:nrow(results)){
  cpg<-rownames(results)[ii]
  weights<-NULL
  betas<-NULL
  studies<-unlist(strsplit(results[ii, "Studies"], ", "))
  for(jj in 1:length(studies)){
    temp<-get(paste(studies[jj]))
    weights<-append(weights, temp[cpg, "weight"])
    betas<-append(betas, temp[cpg, "PTSD"])
    results[cpg, paste(studies[jj], ".beta", sep="")]<-temp[cpg, "PTSD"]
    results[cpg, paste(studies[jj], ".se", sep="")]<-temp[cpg, "s.e."]
    results[cpg, paste(studies[jj], ".N", sep="")]<-temp[cpg, "N.subjects"]
  }
  results[cpg, "beta"]<-sum(betas*weights)/sum(weights)
  results[cpg, "variance"]<-1/sum(weights)
}

data(probe.features)
probe.features<-probe.features[rownames(results), c("CHR", "MAPINFO", "gene", "feature")]
probe.features<-cbind(results, probe.features)
probe.features<-probe.features[, c("Studies", "CpG", "CHR", "MAPINFO", "gene", "feature",
                                   "beta", "variance", "p", "FDR",
                                   paste(studies, ".beta", sep=""),
                                   paste(studies, ".se", sep=""),
                                   paste(studies, ".N", sep=""))]
colnames(probe.features)<-c("Studies", "CpG", "CHR", "Position", "Gene", "Feature",
                            "beta", "variance", "p-value", "FDR",
                            paste(studies, ".beta", sep=""),
                            paste(studies, ".se", sep=""),
                            paste(studies, ".N", sep=""))

rownames(probe.features)<-c(1:nrow(probe.features))
write.csv(probe.features, "PGC_EWAS_wINTRuST_TopResults.csv")

rm(list=ls())

results<-read.csv("PGC_EWAS_wINTRuST_TopResults.csv", row.names=1, stringsAsFactors=F)
studies<-unlist(strsplit(results[1, "Studies"], ","))
studies<-gsub(" ", "", studies)

names<-matrix(nrow=length(studies), ncol=2)
colnames(names)<-c("Old", "Correct")
names[, "Old"]<-studies
rownames(names)<-names[, "Old"]

names[match("VA", rownames(names)), "Correct"]<-"VA-NCPTSD"
names[match("MIR", rownames(names)), "Correct"]<-"VA-M-AA"
names[match("DUKE", rownames(names)), "Correct"]<-"VA-M-EU"
names[match("TRUST", rownames(names)), "Correct"]<-"INTRuST"
names[match("AS", rownames(names)), "Correct"]<-"Army STARRS"
names[which(is.na(names[, "Correct"])), "Correct"]<-names[which(is.na(names[, "Correct"])), "Old"]

pdf("PGC_EWAS_wINTRuST_forestPlots.pdf")
for(ii in 1:nrow(results)){
  betas<-as.numeric(results[ii,paste(studies, ".beta", sep="")])
  stderr<-as.numeric(results[ii,paste(studies, ".se", sep="")])
  lower<-betas-1.96*stderr
  upper<-betas+1.96*stderr
  betaC<-results[ii, "beta"]
  stderrC<-results[ii, "variance"]
  lowerC<-betaC-(1.96*sqrt(stderrC))
  upperC<-betaC+(1.96*sqrt(stderrC))
  beta.table<-c("", "Beta", round(betas,3), NA, round(betaC,3))
  subjs<-results[ii,paste(studies, ".N", sep="")]
  Ns<-as.character(unlist(c("", "N", subjs, NA, sum(subjs))))
  betas<-c(NA, NA, round(betas,3), NA, round(betaC,3))
  lower<-c(NA, NA, lower, NA, lowerC)
  upper<-c(NA, NA, upper, NA, upperC)
  stud<-c(results[ii, "CpG"], "Study", studies, NA, "Summary")
  summary<-cbind(stud, beta.table, Ns)
  summary[which(is.na(summary))]<-""
  # Call forestplot
  cochrane<-data.frame(betas, lower, upper)
  forestplot(summary, cochrane,
             new_page = TRUE,
             is.summary=c(TRUE,TRUE,rep(FALSE,10), TRUE),
             xlog=FALSE, #lineheight=unit(1,"cm"),
             col=fpColors(box="black",line="black", summary="royalblue"),
             xticks=(c(-0.4, -0.2, 0.1)),
             txt_gp = fpTxtGp( label = list(gpar(fontfamily="", cex=1)),
                               ticks = gpar(fontfamily = "", cex=1)))
}
dev.off()

# Significant Sites
sigs<-length(results[as.numeric(results[, "FDR"]<=0.05), "FDR"])
for(ii in 1:sigs){
  betas<-as.numeric(results[ii,paste(studies, ".beta", sep="")])
  stderr<-as.numeric(results[ii,paste(studies, ".se", sep="")])
  lower<-betas-1.96*stderr
  upper<-betas+1.96*stderr
  betaC<-results[ii, "beta"]
  stderrC<-results[ii, "variance"]
  lowerC<-betaC-(1.96*sqrt(stderrC))
  upperC<-betaC+(1.96*sqrt(stderrC))
  beta.table<-c("", "Beta", round(betas,3), NA, round(betaC,3))
  subjs<-results[ii,paste(studies, ".N", sep="")]
  Ns<-as.character(unlist(c("", "N", subjs, NA, sum(subjs))))
  betas<-c(NA, NA, round(betas,3), NA, round(betaC,3))
  lower<-c(NA, NA, lower, NA, lowerC)
  upper<-c(NA, NA, upper, NA, upperC)
  stud<-c(results[ii, "CpG"], "Study", studies, NA, "Summary")
  summary<-cbind(stud, beta.table, Ns)
  summary[which(is.na(summary))]<-""
  
  # Call forestplot
  cochrane<-data.frame(betas, lower, upper)
  png(paste("PGC_EWAS_wINTRuST_", results[ii, "CpG"], ".png", sep=""), width=960, height=640,units="px")
  forestplot(summary, cochrane,
             new_page = TRUE,
             is.summary=c(TRUE,TRUE,rep(FALSE,10),TRUE),
             xlog=FALSE, #lineheight=unit(1,"cm"),
             col=fpColors(box="black",line="black", summary="royalblue"),
             xticks=(c(-0.4, -0.2, 0.1)),
             txt_gp = fpTxtGp( label = list(gpar(fontfamily="", cex=3)),
                               ticks = gpar(fontfamily = "", cex=2)))
  dev.off()
}

# Significant Sites
sigs<-length(results[as.numeric(results[, "FDR"]<=0.05), "FDR"])
names<-names[order(names[, "Correct"]), ]
studies<-names[, "Old"]

for(ii in 1:sigs){
  betas<-as.numeric(results[ii,paste(studies, ".beta", sep="")])
  stderr<-as.numeric(results[ii,paste(studies, ".se", sep="")])
  lower<-betas-1.96*stderr
  upper<-betas+1.96*stderr
  betaC<-results[ii, "beta"]
  stderrC<-results[ii, "variance"]
  lowerC<-betaC-(1.96*sqrt(stderrC))
  upperC<-betaC+(1.96*sqrt(stderrC))
  beta.table<-c("", "Beta", round(betas,3), NA, round(betaC,3))
  subjs<-results[ii,paste(studies, ".N", sep="")]
  Ns<-as.character(unlist(c("", "N", subjs, NA, sum(subjs))))
  betas<-c(NA, NA, round(betas,3), NA, round(betaC,3))
  lower<-c(NA, NA, lower, NA, lowerC)
  upper<-c(NA, NA, upper, NA, upperC)
  stud<-c(results[ii, "CpG"], "Study", names[, "Correct"], NA, "Summary")
  summary<-cbind(stud, beta.table, Ns)
  summary[which(is.na(summary))]<-""
  summary<-summary[, 1]
  
  # Call forestplot
  cochrane<-data.frame(betas, lower, upper)
  summary<-summary[-1]
  cochrane<-cochrane[-1,]
  postscript(paste("PGC_EWAS_wINTRuST_", results[ii, "CpG"], ".eps", sep=""), horizontal=FALSE, width = 3, height = 3, 
             onefile = FALSE, paper = "special", colormodel = "cmyk")
  forestplot(summary, cochrane,
             new_page = TRUE,
             is.summary=c(TRUE, rep(FALSE,10),TRUE),
             xlog=FALSE, #lineheight=unit(1,"cm"),
             col=fpColors(box="grey48",line="grey48", summary="black"),
             xticks=(c(-0.4, -0.2, 0.1)),
             txt_gp = fpTxtGp( label = list(gpar(fontfamily="", cex=0.5)),
                               ticks = gpar(fontfamily = "", cex=0.5)))
  dev.off()

}

# AHRR Sites

names<-names[order(names[, "Correct"]), ]
studies<-names[, "Old"]
cpgs<-c("cg05575921", "cg21161138", "cg25648203", "cg26703534")
results<-results[results[, "CpG"]%in%cpgs, ]

postscript("PGC_EWAS_wINTRuST_AHRR.eps", horizontal=FALSE, width = 6, height = 6, 
           onefile = FALSE, paper = "special", colormodel = "cmyk")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))

for(ii in 1:nrow(results)){
  betas<-as.numeric(results[ii,paste(studies, ".beta", sep="")])
  stderr<-as.numeric(results[ii,paste(studies, ".se", sep="")])
  lower<-betas-1.96*stderr
  upper<-betas+1.96*stderr
  betaC<-results[ii, "beta"]
  stderrC<-results[ii, "variance"]
  lowerC<-betaC-(1.96*sqrt(stderrC))
  upperC<-betaC+(1.96*sqrt(stderrC))
  beta.table<-c("", "Beta", round(betas,3), NA, round(betaC,3))
  subjs<-results[ii,paste(studies, ".N", sep="")]
  Ns<-as.character(unlist(c("", "N", subjs, NA, sum(subjs))))
  betas<-c(NA, NA, round(betas,3), NA, round(betaC,3))
  lower<-c(NA, NA, lower, NA, lowerC)
  upper<-c(NA, NA, upper, NA, upperC)
  stud<-c(results[ii, "CpG"], "Study", names[, "Correct"], NA, "Summary")
  summary<-cbind(stud, beta.table, Ns)
  summary[which(is.na(summary))]<-""
  summary<-summary[, 1]
  
  # Call forestplot
  cochrane<-data.frame(betas, lower, upper)
  summary<-summary[-1]
  cochrane<-cochrane[-1,]
  
  if(ii==1){
    row=1
    col=1
    head(results)
    title=paste("A) ", results[ii, "CpG"], sep="")
  }else if(ii==2){
    row=1
    col=2
    title=paste("B) ", results[ii, "CpG"], sep="")
  }else if(ii==3){
    row=2
    col=1
    title=paste("C) ", results[ii, "CpG"], sep="")
  }else if(ii==4){
    row=2
    col=2
    title=paste("D) ", results[ii, "CpG"], sep="")
  }
  
  pushViewport(viewport(layout.pos.col = col, layout.pos.row = row))
  forestplot(summary, cochrane,
             new_page = FALSE,
             is.summary=c(TRUE, rep(FALSE,10),TRUE),
             xlog=FALSE, #lineheight=unit(1,"cm"),
             col=fpColors(box="grey48",line="grey48", summary="black"),
             xticks=(c(-0.4, -0.2, 0.1)),
             title=title,
             txt_gp = fpTxtGp( label = list(gpar(fontfamily="", cex=0.5)),
                               ticks = gpar(fontfamily = "", cex=0.5)))
  popViewport()
}
dev.off()

png("PGC_EWAS_AHRR.png", width=6, height=6,units="in", res=300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))

for(ii in 1:nrow(results)){
  betas<-as.numeric(results[ii,paste(studies, ".beta", sep="")])
  stderr<-as.numeric(results[ii,paste(studies, ".se", sep="")])
  lower<-betas-1.96*stderr
  upper<-betas+1.96*stderr
  betaC<-results[ii, "beta"]
  stderrC<-results[ii, "variance"]
  lowerC<-betaC-(1.96*sqrt(stderrC))
  upperC<-betaC+(1.96*sqrt(stderrC))
  beta.table<-c("", "Beta", round(betas,3), NA, round(betaC,3))
  subjs<-results[ii,paste(studies, ".N", sep="")]
  Ns<-as.character(unlist(c("", "N", subjs, NA, sum(subjs))))
  betas<-c(NA, NA, round(betas,3), NA, round(betaC,3))
  lower<-c(NA, NA, lower, NA, lowerC)
  upper<-c(NA, NA, upper, NA, upperC)
  stud<-c(results[ii, "CpG"], "Study", names[, "Correct"], NA, "Summary")
  summary<-cbind(stud, beta.table, Ns)
  summary[which(is.na(summary))]<-""
  summary<-summary[, 1]
  
  # Call forestplot
  cochrane<-data.frame(betas, lower, upper)
  summary<-summary[-1]
  cochrane<-cochrane[-1,]

  if(ii==1){
    row=1
    col=1
    head(results)
    title=paste("A) ", results[ii, "CpG"], sep="")
  }else if(ii==2){
    row=1
    col=2
    title=paste("B) ", results[ii, "CpG"], sep="")
  }else if(ii==3){
    row=2
    col=1
    title=paste("C) ", results[ii, "CpG"], sep="")
  }else if(ii==4){
    row=2
    col=2
    title=paste("D) ", results[ii, "CpG"], sep="")
  }
  
  pushViewport(viewport(layout.pos.col = col, layout.pos.row = row))
  forestplot(summary, cochrane,
             new_page = FALSE,
             is.summary=c(TRUE, rep(FALSE,10),TRUE),
             xlog=FALSE, #lineheight=unit(1,"cm"),
             col=fpColors(box="grey48",line="grey48", summary="black"),
             xticks=(c(-0.4, -0.2, 0.1)),
             title=title,
             txt_gp = fpTxtGp( label = list(gpar(fontfamily="", cex=0.5)),
                               ticks = gpar(fontfamily = "", cex=0.5)))
  popViewport()
}
dev.off()
