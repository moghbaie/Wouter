# Mehrnoosh Oghbaie
## HEK293T DMSO/Camptothecin
## U2OS WT/ATM
## 09/21/2020

## Define class and run the project from this file
rm(list=ls())

if (!requireNamespace("rstudioapi", quietly = TRUE))
  install.packages("rstudioapi")

setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

set.seed(123)
CRAN.packages <- c("rstatix","readr","readxl", "xlsx","data.table","reshape2","dplyr","magrittr","fpc","cluster","dendextend","edgeR","gplots","DESeq2",
                   "igraph","sqldf","stringr","corrplot","ggplot2","R6","ggridges","progress","heatmaply","phylogram","Biostrings","ape",
                   "gridExtra","ggrepel","rgl","venn", "writexl","outliers","ggforce", "gtable","factoextra","foreach",'doParallel','BayesFactor')
bioconductor.packages <- c("PTXQC","org.Hs.eg.db","clusterProfiler","ReactomePA", "AnnotationHub","DOSE")

source("functions/Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)
source("functions/Anova_Analysis.R")
source("functions/DE_Analysis.R")
input.dir <- "C:/MQ_projects/2.MQ_runs/Wouter_project/txt"
#runQC(input.dir)


Comparison <- rbind(cbind(rep("HEK293T",4),c(rep(c("CPT",'DMSO'),2)),c(rep("insol",2),rep("wcl",2))),
                    cbind(rep("U2OS",4),c(rep(c("ATM",'WT'),2)),c(rep("insol",2),rep("wcl",2))))

RNAInput <- TemplateRna$
  new()$
  importInput()$
  makeEdgeList()$
  diffExpression()

MSInput <- TemplateProtein$new()$
  importInput(input.dir)$
  removeContaminant()$
  transformData()$
  choosingConditions(Comparison)$
  visualize(Comparison)$
  anovaAnalysis(Comparison)$
  compareRnaProtein(RNAInput)$
  writeFiles()$
  drawScatterplot(RNAInput)


save(MSInput,RNAInput, file='backUp.RData')
