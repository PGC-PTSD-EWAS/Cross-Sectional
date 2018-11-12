########################################################################################
# Limma Analysis for PGC Methylation Meta-Analysis
########################################################################################

rm(list=ls())
library(limma)

### This script runs the following analyses:

# 1. Main model: Dnam = age + gender + Cell Types + Ancestry PCs + PTSD
# 2. Smoking model: Dnam = age + gender + Current Smoking + Cell Types + Ancestry PCs + PTSD
# 3. Main model stratified by gender (if applicable)
# 4. Smoking model stratified by gender (if applicable)
# 5. Main model stratified by smoking status
# 6. Main model stratified by gender and smoking status (if applicable)

###

### Set Working Directory
studyName<-""
setwd("")
if(!dir.exists(paste(getwd(), "/PGC_EWAS_Results/", sep=""))){
  dir.create(paste(getwd(), "/PGC_EWAS_Results/", sep=""))
}
###

########################################################################################
# Step 1: Load  Data and Define Variables
########################################################################################

# Load Phenotype Data
pheno<-read.csv("", row.names=1, stringsAsFactors = F) # Subjects x Variables

# Load Methylation Data
load("") # cpg X subjects
beta.norm<-  # rename DNAm data to beta.norm
# Note: pheno file rows should match beta matrix columns

### Check that the Subjects in the pheno rows matches order of subjects in Beta columns
all(rownames(pheno)==colnames(beta.norm)) # this should be TRUE
###

### Indicate variable names here:

# Notes: 
# if applicable,"gender" should be a character vector with possible responses of "Male" 
# and "Female". The levels should be 1 = "Male" and 2 = "Female". If study is 1 gender 
# only, mark variable as "NA"

# "smoke" should be a factor variable with possible responses of "No" if respondent 
# not a current smoker and "Yes" if respondent is a current a smoker. The levels should
# be 1 = "No" and 2 = "Yes".

# "ptsd" should be a dichotomous factor variable where respondents with a current 
# diagnosis ofPTSD have a value of "Case" and those without a current diagnosis a 
# value of "Control". The levels should be 1 = "Control" and 2 = "Case"

age<-"age"
gender<-"female" # DO NOT DELETE THIS ROW. if study all 1 gender, please put: gender<-NA
PCs<-c("Comp.2", "Comp.3", "Comp.4")
cellTypes<-c("CD8T", "CD4T", "NK", "Bcell", "Mono")
smoke<-"CurrSmoke"
ptsd<-"PTSDpm"

###

# Define Models

if(is.na(gender)){
  main<-as.formula(paste("~", paste(age, paste(PCs, collapse="+"), paste(cellTypes, collapse="+"), ptsd, sep="+"), sep=""))
  smoking<-as.formula(paste("~", paste(age, paste(PCs, collapse="+"), paste(cellTypes, collapse="+"), smoke, ptsd, sep="+"), sep=""))
  mainGenderSpec<-as.formula(paste("~", paste(age, paste(PCs, collapse="+"), paste(cellTypes, collapse="+"), ptsd, sep="+"), sep=""))
  smokingGenderSpec<-as.formula(paste("~", paste(age, paste(PCs, collapse="+"), paste(cellTypes, collapse="+"), smoke, ptsd, sep="+"), sep=""))
  smokingAndGenderSpec<-as.formula(paste("~", paste(age, paste(PCs, collapse="+"), paste(cellTypes, collapse="+"), ptsd, sep="+"), sep=""))
} else{
  main<-as.formula(paste("~", paste(age, gender, paste(PCs, collapse="+"), paste(cellTypes, collapse="+"), ptsd, sep="+"), sep=""))
  smoking<-as.formula(paste("~", paste(age, gender, paste(PCs, collapse="+"), paste(cellTypes, collapse="+"), smoke, ptsd, sep="+"), sep=""))
  mainGenderSpec<-as.formula(paste("~", paste(age, paste(PCs, collapse="+"), paste(cellTypes, collapse="+"), ptsd, sep="+"), sep=""))
  smokingGenderSpec<-as.formula(paste("~", paste(age, paste(PCs, collapse="+"), paste(cellTypes, collapse="+"), smoke, ptsd, sep="+"), sep=""))
  smokingAndGenderSpec<-as.formula(paste("~", paste(age, paste(PCs, collapse="+"), paste(cellTypes, collapse="+"),ptsd, sep="+"), sep=""))
}

models<-c("main", "smoking", "mainGenderSpec", "smokingGenderSpec", "smokingAndGenderSpec")

# Beta values of 0 or 1 will be transformed into M-values of "-Inf". 
range(beta.norm, na.rm=T)

# Changing beta values of 0 to 0.0001
if(min(beta.norm, na.rm=T)==0){
  beta.norm[which(beta.norm==0)]<-0.0001
}

# Changing beta values of 1 to 0.9999
if(max(beta.norm, na.rm=T)==1){
  beta.norm[which(beta.norm==1)]<-0.9999
}

range(beta.norm, na.rm=T)
sum(is.na(beta.norm))

# Convert to Mvalues using log2
beta.norm<-log2(beta.norm/(1-beta.norm)) # log transforming
sum(is.na(beta.norm)) # should be same as the number of missing as above

keep<-c("keep", ls())

########################################################################################
# Step 2: Run Models
########################################################################################

# Main Model
design<-model.matrix(main, data=pheno)
beta<-beta.norm[, rownames(design)]
all(rownames(design)==colnames(beta)) # Should be TRUE
fit<-lmFit(beta, design) # Runs linear models
fit.coef<-fit$coef # Extracts the beta coefficients
all(rownames(fit.coef)==rownames(beta)) # Should be TRUE

p<-pheno[colnames(beta), ]
cases<-rownames(p)[which(p[, ptsd]=="Case")]
controls<-rownames(p)[which(p[, ptsd]=="Control")]
N.subjects<-apply(beta, 1, function(x) sum(!is.na(x)))
Ncases<-apply(beta[, cases], 1, function(x) sum(!is.na(x)))
Ncontrols<-apply(beta[, controls], 1, function(x) sum(!is.na(x)))
meanCases<-apply(beta[, cases], 1, mean, na.rm=T)
meanControls<-apply(beta[, controls], 1, mean, na.rm=T)
sdCases<-apply(beta[, cases], 1, sd, na.rm=T)
sdControls<-apply(beta[, controls], 1, sd, na.rm=T)
fit.coef<-cbind(fit.coef, N.subjects, meanCases, meanControls, sdCases, 
                sdControls, Ncases, Ncontrols)

contrast.matrix<-makeContrasts(paste(ptsd, "Case", sep=""), levels=design)
rownames(contrast.matrix)[1]<-"(Intercept)" # Have to rename the contrast matrix
fit2<-contrasts.fit(fit, contrast.matrix)
fit2.ebayes<-eBayes(fit2) # Run empirical bayes
results<-topTable(fit2.ebayes,number=nrow(fit2.ebayes),coef=1, adjust="BH")
save(fit2.ebayes, fit.coef, results, file=paste(getwd(), "/PGC_EWAS_Results/", studyName, "_main_results.Rdata", sep=""))
rm(list=ls()[-match(keep, ls())])

# Smoking Model
design<-model.matrix(smoking, data=pheno)
beta<-beta.norm[, rownames(design)]
all(rownames(design)==colnames(beta)) # Should be TRUE
fit<-lmFit(beta, design) # Runs linear models
fit.coef<-fit$coef # Extracts the beta coefficients
all(rownames(fit.coef)==rownames(beta)) # Should be TRUE

p<-pheno[colnames(beta), ]
cases<-rownames(p)[which(p[, ptsd]=="Case")]
controls<-rownames(p)[which(p[, ptsd]=="Control")]
N.subjects<-apply(beta, 1, function(x) sum(!is.na(x)))
Ncases<-apply(beta[, cases], 1, function(x) sum(!is.na(x)))
Ncontrols<-apply(beta[, controls], 1, function(x) sum(!is.na(x)))
meanCases<-apply(beta[, cases], 1, mean, na.rm=T)
meanControls<-apply(beta[, controls], 1, mean, na.rm=T)
sdCases<-apply(beta[, cases], 1, sd, na.rm=T)
sdControls<-apply(beta[, controls], 1, sd, na.rm=T)
fit.coef<-cbind(fit.coef, N.subjects, meanCases, meanControls, sdCases, 
                sdControls, Ncases, Ncontrols) 

contrast.matrix<-makeContrasts(paste(ptsd, "Case", sep=""), levels=design)
rownames(contrast.matrix)[1]<-"(Intercept)" # Have to rename the contrast matrix
fit2<-contrasts.fit(fit, contrast.matrix)
fit2.ebayes<-eBayes(fit2) # Run empirical bayes
results<-topTable(fit2.ebayes,number=nrow(fit2.ebayes),coef=1, adjust="BH")
save(fit2.ebayes, fit.coef, results, file=paste(getwd(), "/PGC_EWAS_Results/", studyName, "_smoking_results.Rdata", sep=""))
rm(list=ls()[-match(keep, ls())])

# Main Gender Specific
if(!is.na(gender)){
  # Males
  p<-pheno[pheno[, gender]=="Male", ]
  design<-model.matrix(mainGenderSpec, data=p)
  rm(p)
  beta<-beta.norm[, rownames(design)]
  all(rownames(design)==colnames(beta)) # Should be TRUE
  fit<-lmFit(beta, design) # Runs linear models
  fit.coef<-fit$coef # Extracts the beta coefficients
  all(rownames(fit.coef)==rownames(beta)) # Should be TRUE
  
  p<-pheno[colnames(beta), ]
  cases<-rownames(p)[which(p[, ptsd]=="Case")]
  controls<-rownames(p)[which(p[, ptsd]=="Control")]
  N.subjects<-apply(beta, 1, function(x) sum(!is.na(x)))
  Ncases<-apply(beta[, cases], 1, function(x) sum(!is.na(x)))
  Ncontrols<-apply(beta[, controls], 1, function(x) sum(!is.na(x)))
  meanCases<-apply(beta[, cases], 1, mean, na.rm=T)
  meanControls<-apply(beta[, controls], 1, mean, na.rm=T)
  sdCases<-apply(beta[, cases], 1, sd, na.rm=T)
  sdControls<-apply(beta[, controls], 1, sd, na.rm=T)
  fit.coef<-cbind(fit.coef, N.subjects, meanCases, meanControls, sdCases, 
                  sdControls, Ncases, Ncontrols)
  
  contrast.matrix<-makeContrasts(paste(ptsd, "Case", sep=""), levels=design)
  rownames(contrast.matrix)[1]<-"(Intercept)" # Have to rename the contrast matrix
  fit2<-contrasts.fit(fit, contrast.matrix)
  fit2.ebayes<-eBayes(fit2) # Run empirical bayes
  results<-topTable(fit2.ebayes,number=nrow(fit2.ebayes),coef=1, adjust="BH")
  save(fit2.ebayes, fit.coef, results, file=paste(getwd(), "/PGC_EWAS_Results/", studyName, "_main_results_malesOnly.Rdata", sep=""))
  rm(list=ls()[-match(keep, ls())])
  
  # Females
  p<-pheno[pheno[, gender]=="Female", ]
  design<-model.matrix(mainGenderSpec, data=p)
  rm(p)
  beta<-beta.norm[, rownames(design)]
  all(rownames(design)==colnames(beta)) # Should be TRUE
  fit<-lmFit(beta, design) # Runs linear models
  fit.coef<-fit$coef # Extracts the beta coefficients
  all(rownames(fit.coef)==rownames(beta)) # Should be TRUE
  
  p<-pheno[colnames(beta), ]
  cases<-rownames(p)[which(p[, ptsd]=="Case")]
  controls<-rownames(p)[which(p[, ptsd]=="Control")]
  N.subjects<-apply(beta, 1, function(x) sum(!is.na(x)))
  Ncases<-apply(beta[, cases], 1, function(x) sum(!is.na(x)))
  Ncontrols<-apply(beta[, controls], 1, function(x) sum(!is.na(x)))
  meanCases<-apply(beta[, cases], 1, mean, na.rm=T)
  meanControls<-apply(beta[, controls], 1, mean, na.rm=T)
  sdCases<-apply(beta[, cases], 1, sd, na.rm=T)
  sdControls<-apply(beta[, controls], 1, sd, na.rm=T)
  fit.coef<-cbind(fit.coef, N.subjects, meanCases, meanControls, sdCases, 
                  sdControls, Ncases, Ncontrols)
  
  contrast.matrix<-makeContrasts(paste(ptsd, "Case", sep=""), levels=design)
  rownames(contrast.matrix)[1]<-"(Intercept)" # Have to rename the contrast matrix
  fit2<-contrasts.fit(fit, contrast.matrix)
  fit2.ebayes<-eBayes(fit2) # Run empirical bayes
  results<-topTable(fit2.ebayes,number=nrow(fit2.ebayes),coef=1, adjust="BH")
  save(fit2.ebayes, fit.coef, results, file=paste(getwd(), "/PGC_EWAS_Results/", studyName, "_main_results_femalesOnly.Rdata", sep=""))
  rm(list=ls()[-match(keep, ls())])
}

# Smoking Gender Specific
if(!is.na(gender)){
  # Males
  p<-pheno[pheno[, gender]=="Male", ]
  design<-model.matrix(smokingGenderSpec, data=p)
  rm(p)
  beta<-beta.norm[, rownames(design)]
  all(rownames(design)==colnames(beta)) # Should be TRUE
  fit<-lmFit(beta, design) # Runs linear models
  fit.coef<-fit$coef # Extracts the beta coefficients
  all(rownames(fit.coef)==rownames(beta)) # Should be TRUE
  
  p<-pheno[colnames(beta), ]
  cases<-rownames(p)[which(p[, ptsd]=="Case")]
  controls<-rownames(p)[which(p[, ptsd]=="Control")]
  N.subjects<-apply(beta, 1, function(x) sum(!is.na(x)))
  Ncases<-apply(beta[, cases], 1, function(x) sum(!is.na(x)))
  Ncontrols<-apply(beta[, controls], 1, function(x) sum(!is.na(x)))
  meanCases<-apply(beta[, cases], 1, mean, na.rm=T)
  meanControls<-apply(beta[, controls], 1, mean, na.rm=T)
  sdCases<-apply(beta[, cases], 1, sd, na.rm=T)
  sdControls<-apply(beta[, controls], 1, sd, na.rm=T)
  fit.coef<-cbind(fit.coef, N.subjects, meanCases, meanControls, sdCases, 
                  sdControls, Ncases, Ncontrols)
  
  contrast.matrix<-makeContrasts(paste(ptsd, "Case", sep=""), levels=design)
  rownames(contrast.matrix)[1]<-"(Intercept)" # Have to rename the contrast matrix
  fit2<-contrasts.fit(fit, contrast.matrix)
  fit2.ebayes<-eBayes(fit2) # Run empirical bayes
  results<-topTable(fit2.ebayes,number=nrow(fit2.ebayes),coef=1, adjust="BH")
  save(fit2.ebayes, fit.coef, results, file=paste(getwd(), "/PGC_EWAS_Results/", studyName, "_smoking_results_malesOnly.Rdata", sep=""))
  rm(list=ls()[-match(keep, ls())])
  
  # Females
  p<-pheno[pheno[, gender]=="Female", ]
  design<-model.matrix(smokingGenderSpec, data=p)
  
  beta<-beta.norm[, rownames(design)]
  all(rownames(design)==colnames(beta)) # Should be TRUE
  fit<-lmFit(beta, design) # Runs linear models
  fit.coef<-fit$coef # Extracts the beta coefficients
  all(rownames(fit.coef)==rownames(beta)) # Should be TRUE
  
  p<-pheno[colnames(beta), ]
  cases<-rownames(p)[which(p[, ptsd]=="Case")]
  controls<-rownames(p)[which(p[, ptsd]=="Control")]
  N.subjects<-apply(beta, 1, function(x) sum(!is.na(x)))
  Ncases<-apply(beta[, cases], 1, function(x) sum(!is.na(x)))
  Ncontrols<-apply(beta[, controls], 1, function(x) sum(!is.na(x)))
  meanCases<-apply(beta[, cases], 1, mean, na.rm=T)
  meanControls<-apply(beta[, controls], 1, mean, na.rm=T)
  sdCases<-apply(beta[, cases], 1, sd, na.rm=T)
  sdControls<-apply(beta[, controls], 1, sd, na.rm=T)
  fit.coef<-cbind(fit.coef, N.subjects, meanCases, meanControls, sdCases, 
                  sdControls, Ncases, Ncontrols)
  
  contrast.matrix<-makeContrasts(paste(ptsd, "Case", sep=""), levels=design)
  rownames(contrast.matrix)[1]<-"(Intercept)" # Have to rename the contrast matrix
  fit2<-contrasts.fit(fit, contrast.matrix)
  fit2.ebayes<-eBayes(fit2) # Run empirical bayes
  results<-topTable(fit2.ebayes,number=nrow(fit2.ebayes),coef=1, adjust="BH")
  save(fit2.ebayes, fit.coef, results, file=paste(getwd(), "/PGC_EWAS_Results/", studyName, "_smoking_results_femalesOnly.Rdata", sep=""))
  rm(list=ls()[-match(keep, ls())])
}

# Main Model among Non-Smokers
p<-pheno[which(pheno[, smoke]=="No"),]
design<-model.matrix(main, data=p)
rm(p)
beta<-beta.norm[, rownames(design)]
all(rownames(design)==colnames(beta)) # Should be TRUE
fit<-lmFit(beta, design) # Runs linear models
fit.coef<-fit$coef # Extracts the beta coefficients
all(rownames(fit.coef)==rownames(beta)) # Should be TRUE

p<-pheno[colnames(beta), ]
cases<-rownames(p)[which(p[, ptsd]=="Case")]
controls<-rownames(p)[which(p[, ptsd]=="Control")]
N.subjects<-apply(beta, 1, function(x) sum(!is.na(x)))
Ncases<-apply(beta[, cases], 1, function(x) sum(!is.na(x)))
Ncontrols<-apply(beta[, controls], 1, function(x) sum(!is.na(x)))
meanCases<-apply(beta[, cases], 1, mean, na.rm=T)
meanControls<-apply(beta[, controls], 1, mean, na.rm=T)
sdCases<-apply(beta[, cases], 1, sd, na.rm=T)
sdControls<-apply(beta[, controls], 1, sd, na.rm=T)
fit.coef<-cbind(fit.coef, N.subjects, meanCases, meanControls, sdCases, 
                sdControls, Ncases, Ncontrols)

contrast.matrix<-makeContrasts(paste(ptsd, "Case", sep=""), levels=design)
rownames(contrast.matrix)[1]<-"(Intercept)" # Have to rename the contrast matrix
fit2<-contrasts.fit(fit, contrast.matrix)
fit2.ebayes<-eBayes(fit2) # Run empirical bayes
results<-topTable(fit2.ebayes,number=nrow(fit2.ebayes),coef=1, adjust="BH")
save(fit2.ebayes, fit.coef, results, file=paste(getwd(), "/PGC_EWAS_Results/", studyName, "_main_results_nonSmokersOnly.Rdata", sep=""))
rm(list=ls()[-match(keep, ls())])

# Main Model among Smokers Only
p<-pheno[which(pheno[, smoke]=="Yes"),]
design<-model.matrix(mainGenderSpec, data=p)
rm(p)
beta<-beta.norm[, rownames(design)]
all(rownames(design)==colnames(beta)) # Should be TRUE
fit<-lmFit(beta, design) # Runs linear models
fit.coef<-fit$coef # Extracts the beta coefficients
all(rownames(fit.coef)==rownames(beta)) # Should be TRUE

p<-pheno[colnames(beta), ]
cases<-rownames(p)[which(p[, ptsd]=="Case")]
controls<-rownames(p)[which(p[, ptsd]=="Control")]
N.subjects<-apply(beta, 1, function(x) sum(!is.na(x)))
Ncases<-apply(beta[, cases], 1, function(x) sum(!is.na(x)))
Ncontrols<-apply(beta[, controls], 1, function(x) sum(!is.na(x)))
meanCases<-apply(beta[, cases], 1, mean, na.rm=T)
meanControls<-apply(beta[, controls], 1, mean, na.rm=T)
sdCases<-apply(beta[, cases], 1, sd, na.rm=T)
sdControls<-apply(beta[, controls], 1, sd, na.rm=T)
fit.coef<-cbind(fit.coef, N.subjects, meanCases, meanControls, sdCases, 
                sdControls, Ncases, Ncontrols)

contrast.matrix<-makeContrasts(paste(ptsd, "Case", sep=""), levels=design)
rownames(contrast.matrix)[1]<-"(Intercept)" # Have to rename the contrast matrix
fit2<-contrasts.fit(fit, contrast.matrix)
fit2.ebayes<-eBayes(fit2) # Run empirical bayes
results<-topTable(fit2.ebayes,number=nrow(fit2.ebayes),coef=1, adjust="BH")
save(fit2.ebayes, fit.coef, results, file=paste(getwd(), "/PGC_EWAS_Results/", studyName, "_main_results_smokersOnly.Rdata", sep=""))
rm(list=ls()[-match(keep, ls())])



# Main Model among Non-Smokers Gender Specific
if(!is.na(gender)){
  # Non-Smoking Males
  p<-pheno[which(pheno[, gender]=="Male" & pheno[, smoke]=="No"), ]
  design<-model.matrix(mainGenderSpec, data=p)
  rm(p)
  beta<-beta.norm[, rownames(design)]
  all(rownames(design)==colnames(beta)) # Should be TRUE
  fit<-lmFit(beta, design) # Runs linear models
  fit.coef<-fit$coef # Extracts the beta coefficients
  all(rownames(fit.coef)==rownames(beta)) # Should be TRUE
  
  p<-pheno[colnames(beta), ]
  cases<-rownames(p)[which(p[, ptsd]=="Case")]
  controls<-rownames(p)[which(p[, ptsd]=="Control")]
  N.subjects<-apply(beta, 1, function(x) sum(!is.na(x)))
  Ncases<-apply(beta[, cases], 1, function(x) sum(!is.na(x)))
  Ncontrols<-apply(beta[, controls], 1, function(x) sum(!is.na(x)))
  meanCases<-apply(beta[, cases], 1, mean, na.rm=T)
  meanControls<-apply(beta[, controls], 1, mean, na.rm=T)
  sdCases<-apply(beta[, cases], 1, sd, na.rm=T)
  sdControls<-apply(beta[, controls], 1, sd, na.rm=T)
  fit.coef<-cbind(fit.coef, N.subjects, meanCases, meanControls, sdCases, 
                  sdControls, Ncases, Ncontrols)
  
  contrast.matrix<-makeContrasts(paste(ptsd, "Case", sep=""), levels=design)
  rownames(contrast.matrix)[1]<-"(Intercept)" # Have to rename the contrast matrix
  fit2<-contrasts.fit(fit, contrast.matrix)
  fit2.ebayes<-eBayes(fit2) # Run empirical bayes
  results<-topTable(fit2.ebayes,number=nrow(fit2.ebayes),coef=1, adjust="BH")
  save(fit2.ebayes, fit.coef, results, file=paste(getwd(), "/PGC_EWAS_Results/", studyName, "_main_results_maleNonSmokersOnly.Rdata", sep=""))
  rm(list=ls()[-match(keep, ls())])
  
  # Non-Smoking Females
  p<-pheno[which(pheno[, gender]=="Female" & pheno[, smoke]=="No"), ]
  design<-model.matrix(mainGenderSpec, data=p)
  rm(p)
  beta<-beta.norm[, rownames(design)]
  all(rownames(design)==colnames(beta)) # Should be TRUE
  fit<-lmFit(beta, design) # Runs linear models
  fit.coef<-fit$coef # Extracts the beta coefficients
  all(rownames(fit.coef)==rownames(beta)) # Should be TRUE
  
  p<-pheno[colnames(beta), ]
  cases<-rownames(p)[which(p[, ptsd]=="Case")]
  controls<-rownames(p)[which(p[, ptsd]=="Control")]
  N.subjects<-apply(beta, 1, function(x) sum(!is.na(x)))
  Ncases<-apply(beta[, cases], 1, function(x) sum(!is.na(x)))
  Ncontrols<-apply(beta[, controls], 1, function(x) sum(!is.na(x)))
  meanCases<-apply(beta[, cases], 1, mean, na.rm=T)
  meanControls<-apply(beta[, controls], 1, mean, na.rm=T)
  sdCases<-apply(beta[, cases], 1, sd, na.rm=T)
  sdControls<-apply(beta[, controls], 1, sd, na.rm=T)
  fit.coef<-cbind(fit.coef, N.subjects, meanCases, meanControls, sdCases, 
                  sdControls, Ncases, Ncontrols)
  
  contrast.matrix<-makeContrasts(paste(ptsd, "Case", sep=""), levels=design)
  rownames(contrast.matrix)[1]<-"(Intercept)" # Have to rename the contrast matrix
  fit2<-contrasts.fit(fit, contrast.matrix)
  fit2.ebayes<-eBayes(fit2) # Run empirical bayes
  results<-topTable(fit2.ebayes,number=nrow(fit2.ebayes),coef=1, adjust="BH")
  save(fit2.ebayes, fit.coef, results, file=paste(getwd(), "/PGC_EWAS_Results/", studyName, "_main_results_femaleNonSmokersOnly.Rdata", sep=""))
  rm(list=ls()[-match(keep, ls())])
  
  # Smoking Females
  p<-pheno[which(pheno[, gender]=="Female" & pheno[, smoke]=="Yes"), ]
  design<-model.matrix(mainGenderSpec, data=p)
  rm(p)
  beta<-beta.norm[, rownames(design)]
  all(rownames(design)==colnames(beta)) # Should be TRUE
  fit<-lmFit(beta, design) # Runs linear models
  fit.coef<-fit$coef # Extracts the beta coefficients
  all(rownames(fit.coef)==rownames(beta)) # Should be TRUE
  
  p<-pheno[colnames(beta), ]
  cases<-rownames(p)[which(p[, ptsd]=="Case")]
  controls<-rownames(p)[which(p[, ptsd]=="Control")]
  N.subjects<-apply(beta, 1, function(x) sum(!is.na(x)))
  Ncases<-apply(beta[, cases], 1, function(x) sum(!is.na(x)))
  Ncontrols<-apply(beta[, controls], 1, function(x) sum(!is.na(x)))
  meanCases<-apply(beta[, cases], 1, mean, na.rm=T)
  meanControls<-apply(beta[, controls], 1, mean, na.rm=T)
  sdCases<-apply(beta[, cases], 1, sd, na.rm=T)
  sdControls<-apply(beta[, controls], 1, sd, na.rm=T)
  fit.coef<-cbind(fit.coef, N.subjects, meanCases, meanControls, sdCases, 
                  sdControls, Ncases, Ncontrols)
  
  contrast.matrix<-makeContrasts(paste(ptsd, "Case", sep=""), levels=design)
  rownames(contrast.matrix)[1]<-"(Intercept)" # Have to rename the contrast matrix
  fit2<-contrasts.fit(fit, contrast.matrix)
  fit2.ebayes<-eBayes(fit2) # Run empirical bayes
  results<-topTable(fit2.ebayes,number=nrow(fit2.ebayes),coef=1, adjust="BH")
  save(fit2.ebayes, fit.coef, results, file=paste(getwd(), "/PGC_EWAS_Results/", studyName, "_main_results_maleSmokersOnly.Rdata", sep=""))
  rm(list=ls()[-match(keep, ls())])
  
  p<-pheno[which(pheno[, gender]=="Female" & pheno[, smoke]=="No"), ]
  design<-model.matrix(mainGenderSpec, data=p)
  rm(p)
  beta<-beta.norm[, rownames(design)]
  all(rownames(design)==colnames(beta)) # Should be TRUE
  fit<-lmFit(beta, design) # Runs linear models
  fit.coef<-fit$coef # Extracts the beta coefficients
  all(rownames(fit.coef)==rownames(beta)) # Should be TRUE
  
  p<-pheno[colnames(beta), ]
  cases<-rownames(p)[which(p[, ptsd]=="Case")]
  controls<-rownames(p)[which(p[, ptsd]=="Control")]
  N.subjects<-apply(beta, 1, function(x) sum(!is.na(x)))
  Ncases<-apply(beta[, cases], 1, function(x) sum(!is.na(x)))
  Ncontrols<-apply(beta[, controls], 1, function(x) sum(!is.na(x)))
  meanCases<-apply(beta[, cases], 1, mean, na.rm=T)
  meanControls<-apply(beta[, controls], 1, mean, na.rm=T)
  sdCases<-apply(beta[, cases], 1, sd, na.rm=T)
  sdControls<-apply(beta[, controls], 1, sd, na.rm=T)
  fit.coef<-cbind(fit.coef, N.subjects, meanCases, meanControls, sdCases, 
                  sdControls, Ncases, Ncontrols)
  
  contrast.matrix<-makeContrasts(paste(ptsd, "Case", sep=""), levels=design)
  rownames(contrast.matrix)[1]<-"(Intercept)" # Have to rename the contrast matrix
  fit2<-contrasts.fit(fit, contrast.matrix)
  fit2.ebayes<-eBayes(fit2) # Run empirical bayes
  results<-topTable(fit2.ebayes,number=nrow(fit2.ebayes),coef=1, adjust="BH")
  save(fit2.ebayes, fit.coef, results, file=paste(getwd(), "/PGC_EWAS_Results/", studyName, "_main_results_femaleSmokersOnly.Rdata", sep=""))
  rm(list=ls()[-match(keep, ls())])
}




