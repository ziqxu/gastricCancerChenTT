while(!require(edgeR)){install.packages('edgeR')}
while(!require(magrittr)){install.packages('magrittr')}
while(!require(biomaRt)){install.packages('biomaRt')}
while(!require(tibble)){install.packages('tibble')}
while(!require(dplyr)){install.packages('dplyr')}
library(xlsx)
library(tibble)
library(survival)
library(rms)
library(glmnet)
library(randomForest)

{
  library(survminer)
  library(gridExtra)
  library(ggpmisc)
  library(ggpubr)
  library(grid)
  library(rmda)
  library(survivalROC)
}