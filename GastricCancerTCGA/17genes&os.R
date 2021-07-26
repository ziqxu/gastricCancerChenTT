source('load_package.R')
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
  result <- ifelse(x<median_x,0,1)
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

full.anal <- function(){
  y <- coxph(Surv(time = OS.time,event = OS)~.-OS-OS.time,data = os.anal.df)%>%summary()
  y <- cbind(y$coefficient, y$conf.int)
  y <- y[,c(2,8,9,5),drop=F]
  # y<-coefficients(y)
  return(y)
}
mul.result <- full.anal()
mul.result

#univariate cox analysis
# uni.anal <- function(x){
#   fmt.str <- paste0('Surv(time = OS.time,event = OS)~',x)%>%
#     as.formula()
#   y <- coxph(fmt.str, data = os.anal.df)%>%summary()
#   y <- cbind(y$coefficient, y$conf.int)
#   y <- y[1,c(2,8,9,5),drop=F]
#   return(y)
# }
# uni.result<-do.call(rbind,lapply(gene.list,uni.anal))
# result.table <- cbind(uni.result,mul.result)


####rms packages part
fml<-paste0('Surv(time = OS.time,event = OS)~',paste0(gene.list,collapse = '+'))%>%as.formula()
coxm<-cph(fml,data = os.anal.df,x=T,y=T,surv = T)
med<-Quantile(coxm)
surv<-Survival(coxm)
surv1<-function(x) surv(times = 30*12,lp = x)
surv2<-function(x) surv(times = 30*36,lp = x)
surv3<-function(x) surv(times = 30*60,lp = x)
dd <- datadist(os.anal.df)
options(datadist=dd)
nom <- nomogram(coxm,fun=list(surv1,surv2,surv3),lp=F,funlabel = c("1-year OS","3-year OS",'5-year OS'))
plot(nom,xfrac = 0.3)

nomogram(coxm,function(x) med(lp=x))


