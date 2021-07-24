load('gastric_tcga_tmm.RData')
rm(list=ls()[!ls()%in%c('tmm.hugo.df')])
clin.df<-read.table('TCGA-STAD.survival.tsv',header = T)
library(xlsx)
library(tibble)
gene.list<-xlsx::read.xlsx('17 genes(1).xlsx',sheetIndex = 1,header = F)%>%.$X1

if(all(gene.list%in%tmm.hugo.df$hgnc_symbol)){
  gene.list.df <- tmm.hugo.df[which(tmm.hugo.df$hgnc_symbol %in% gene.list),]%>%
    .[,-1]%>%
    remove_rownames()%>%
    column_to_rownames(var = 'hgnc_symbol')%>%
    t()%>%
    as.data.frame()%>%
    rownames_to_column(var = 'sample')
}
rm(tmm.hugo.df)

gene.list.df$sample%<>%gsub('\\.','-',x=.)

os.df<-merge(gene.list.df,clin.df,by = 'sample')%>%
  mutate(typecode=substr(sample,nchar(sample)-2,nchar(sample)))
table(os.df$typecode)

os.df%<>%filter(typecode%in%c('01A','01B'))

###clear for space
rm(list = ls()[!ls()%in%c('os.df','gene.list')])

###median grouping
rep.df<-os.df
groupby_median <- function(x){
  median_x <- median(x)
  result <- ifelse(x>median_x,'high','low')
  return(result)
}
rep.df[gene.list] <- lapply(rep.df[gene.list],groupby_median)
cat.df<-rep.df
rm(rep.df)

if(all(!duplicated(cat.df$X_PATIENT))){
  os.anal.df <- cat.df[c(gene.list,'OS','OS.time')]
}
###load survival analysis packages
library(survival)

#multivariate cox analysis
coxph(Surv(time = OS.time,event = OS)~.-OS-OS.time,data = os.anal.df)


#univariate cox analysis

uni.anal <- function(x){
  fmt.str <- paste0('Surv(time = OS.time,event = OS)~',x)%>%
    as.formula()
  y <- coxph(fmt.str, data = os.anal.df)
  return(y)
}
lapply(gene.list,uni.anal)


