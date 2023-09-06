library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(wateRmelon)
library(future)
library(gtools)
library(matrixStats)
library(data.table)
library(MASS)
library(sandwich)
library(lmtest)
library(parallel)
library(R.utils)
library(data.table)
library(readxl)
library(matrixStats)
library(MASS)
library(dplyr)
library(tidyverse)
library(compare)

beta2 <- readRDS("../../1qc6/2000HIV.Mvalue.rds")

beta3 <- beta2[ ,-c(318,325,328)]

pheno3=readRDS("../../2Pheno/Phe_dis_clean_PCs.rds")
pheno3$ID_idat <- pheno3$ID

covid <- read.csv("covid_groups_selected_variables_may.csv", sep=",")

pheno4 <- merge(pheno3, covid, by = "ID", all=T)
pheno4 <- dplyr::filter(pheno4, !is.na(pheno4$ID_idat))
pheno4 <- dplyr::filter(pheno4, pheno4$ETHNICITY=="White")
pheno5 <- pheno4 %>% dplyr::select("ID", "AGE", "SEX_BIRTH", "Sample_Plate", "season_sin", "season_cos", "cov", "vac", "iso2", "cplv", "SMOKING", "BMI_BASELINE", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu", "PC1", "PC2", "PC3", "PC4", "PC5")

pheno5 <- dplyr::filter(pheno5, !is.na(pheno5$vac))
idx <- intersect(colnames(beta3), pheno5$ID)
beta3 <- beta3[, idx]
idx2 <- match(colnames(beta3), pheno5$ID)
pheno6 <- pheno5[idx2, ]

compare(colnames(beta3), pheno6$ID)
cat(dim(pheno6), dim(beta3))

EWAS_pheno=data.frame()
for (i in 1:length(pheno6$ID)) {
    id = pheno6$ID[i]
    #cat(i, id, "|")
    EWAS_pheno[i, "ID_idat"] = id
    INR = dplyr::filter(pheno6, pheno6$ID==id)[, "vac"] 
    if (is.na(INR)) {
        EWAS_pheno[i, "Control"] = NA
    } else {
        INR = as.numeric(INR)
        EWAS_pheno[i, "Control"] = INR
    }
}

EWAS_pheno$Control <- as.factor(EWAS_pheno$Control)
ix <- (which(is.na(EWAS_pheno$Control)))
if (length(ix)>0) {
    beta3 <- beta3[ ,-ix]
    EWAS_pheno <- EWAS_pheno[-ix, ]
    pheno6 <- pheno6[-ix, ]
}

cat("Adter removing the NA individuals. Do the Pheno5 and EWAS_Pheno have the same ID_idat: ")
compare(pheno6$ID, EWAS_pheno$ID_idat)
cat("Adter removing the NA individuals. Do the beta3's colnames are the same as ID_idat in pheno5 and EWAS_pheno: ")
compare(colnames(beta3), EWAS_pheno$ID_idat)
cat("Adter removing the NA individuals. The dimentions of EWAS_pheno, beta3 and pheno5: ")
cat(dim(EWAS_pheno), dim(beta3), dim(pheno6))

removeOutliers<-function(probes){
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR 
  maskU <- probes > row2575[,2] + 3 * rowIQR 
  initial_NAs<-rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe<-rowSums(!is.na(probes))
  Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  return(list(probes, Log))
}

#Remove outliers from METH (methylation data where probes are rows and samples are columns)
system.time(OutlierResults<-removeOutliers(beta3))  
METH.2<-OutlierResults[[1]]
Log<-OutlierResults[[2]]
rm(beta2)
rm(beta3)
save(Log,file="Outlier_log_vacEur.Rdata") #save log

mVals.T0=t(METH.2)
pheno.T0=pheno6
TrImmRes.T0=EWAS_pheno

library("tidyr")
#https://sparkbyexamples.com/r-programming/remove-rows-with-na-in-r/

genes <- unique(colnames(mVals.T0))#test remove [1:100]
df=mVals.T0 
df = df %>% as.data.frame()
gene=genes

cor_test=function(gene)
{ 
  sub <- df[, gene] %>% data.frame()  #df=normalized gene count matrix
  colnames(sub) <- 'methylation'
  test.df <- cbind(pheno.T0, sub)
  cyto1 <- TrImmRes.T0[,2]
  bad <- as.numeric(rep(NA, 4))
  names(bad) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  result <- bad  
  tryCatch(
      {
        ML = rlm(methylation ~ cyto1 + AGE + SEX_BIRTH + Sample_Plate + season_sin + season_cos +
               CD8T + CD4T + NK + Bcell + Mono + Neu, 
               data=test.df)
        cf <- try(coeftest(ML, vcov=vcovHC(ML, type="HC0")))
        x <<- x+1
        a <- sum(x:length(P_value_Index))
        cat(x, ":",length(P_value_Index),".\n",sep="")
        if (class(cf)=="try-error") {
          bad <- as.numeric(rep(NA, 4))
          names(bad)<- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
          result <- bad
        }
        else{
          result <- cf[2, c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
        }  
      },
      error=function(error_message) {
        message("This is my custom message.")
        message("And below is the error message from R:")
        message(error_message)
        return(result)
      }
    )
  return(result)
}

pheno.T0$SEX_BIRTH = as.factor(pheno.T0$SEX_BIRTH)
pheno.T0$Sample_Plate = droplevels(as.factor(pheno.T0$Sample_Plate))
pheno.T0$AGE = as.numeric(pheno.T0$AGE)

L = list()
P_value_Index=genes
cor_info <- data.frame()
x <- 0
Result=vapply(P_value_Index, cor_test, numeric(4))
cor_info=as.data.frame(t(Result))
colnames(cor_info)=c("Estimate", "StdError", "z-score", "pval")
cor_info[,5]=rownames(cor_info)#cor_info[,4]=genes
colnames(cor_info)=c("Estimate", "StdError", "z-score", "pval", "CpGsite")
#saveRDS(cor_info,"cor_info.EWAS_CANNABIS_FREQUENT_INLIFETIME_zscore_noNA.rds",compress="xz")
cor_info$FDR=p.adjust(cor_info$pval, method = "fdr")
saveRDS(cor_info,"corEWASrlm_dis_vac_Eur.rds",compress="xz")



