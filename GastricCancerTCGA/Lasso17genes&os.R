source('load_package.R')
load('gastric_tcga_tmm.RData')
rm(list=ls()[!ls()%in%c('tmm.hugo.df')])

clin.df<-read.table('TCGA-STAD.survival.tsv',header = T)
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
rownames(os.df)<-os.df$sample
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

#multivariate cox analysis
seed <- 101
set.seed(seed)
# train.set <- sample(1:nrow(os.anal.df),nrow(os.anal.df))
# test.set <- setdiff(1:nrow(os.anal.df),train.set)
train.set <- 1:nrow(os.anal.df)

x.mat <- os.anal.df[train.set,gene.list]%>%as.matrix()
y.mat <- os.anal.df[train.set,c('OS.time','OS')]%>%as.matrix()
surv <- Surv(time = os.anal.df$OS.time[train.set],event = os.anal.df$OS[train.set])

fit <- glmnet(x.mat,surv,family = 'cox')
set.seed(seed)
cv.fit <- cv.glmnet(x.mat,surv,family = 'cox',nfolds = 10)
cv.fit$lambda.min

update.list <- coef(fit,s = cv.fit$lambda.min)%>%as.matrix()%>%as.data.frame%>%filter(`1`!=0)%>%rownames()
#update.list <- c("ADAM12","CDH19","SLC5A7","SHOX2","C1QL1","RHPN2","APLN","MNX1","IGF2BP1","HMGB3")
#update.list <- c("ADAM12","SHOX2","MNX1")
fml<-paste0('surv~',paste0(update.list,collapse = '+'))%>%as.formula()
wt.coef<-coef(coxph(fml,data = os.anal.df[train.set,]))%>%as.matrix()
risk.score<-x.mat[,update.list]%*%wt.coef
median(risk.score)

train.df <- os.anal.df[train.set,]
train.df$score <- risk.score
####################
full.anal <- function(x){
  surv<-Surv(time = os.anal.df$OS.time,event = os.anal.df$OS)
  fml<-paste0('surv~',paste0(x,collapse = '+'))%>%as.formula()
  y <- coxph(fml,data = os.anal.df)%>%summary()
  # y <- cbind(y$coefficient, y$conf.int)
  # y <- y[,c(1,2,8,9,5),drop=F]
  # y<-coefficients(y)
  return(y)
}
mul.result <- full.anal(update.list)
mul.result


#univariate cox analysis
uni.anal <- function(x,df){
  fmt.str <- paste0('Surv(time = OS.time,event = OS)~',x)%>%
    as.formula()
  y <- coxph(fmt.str, data = df)%>%summary()
  y <- cbind(y$coefficient, y$conf.int,concordance=y$concordance[1],c.se=y$concordance[2])
  y <- y[1,c(2,8,9,5,10,11),drop=F]
  
  return(y)
}

uni.anal('score',train.df)
uni.result<-do.call(rbind,lapply(gene.list,function(x)uni.anal(x,train.df)))
# result.table <- cbind(uni.result,mul.result)
train.df$score.group<-ifelse(train.df$score>=median(risk.score),'high',"low")%>%
  as.factor()%>%relevel(ref = 'low')
uni.anal('score.group',train.df)
####Validation group
# test.x.mat <- as.matrix(os.anal.df[test.set,])
# test.df <- os.anal.df[test.set,]
# test.df$score<-test.x.mat[,update.list]%*%wt.coef
# test.df$score.group <- ifelse(train.df$score>=median(test.df$score),'high',"low")%>%
#   as.factor()%>%relevel(ref = 'low')
# table(test.df$score.group)
# uni.anal('score.group',test.df)



####rms packages part

pheno<-readr::read_tsv('TCGA-STAD.GDC_phenotype.tsv')%>%
  filter(submitter_id.samples%in%rownames(train.df))%>%
  column_to_rownames(var='submitter_id.samples')
mulcox.df<-merge(train.df,pheno, by = 0)
mulcox.list <- c('score.group',)

fml<-paste0('Surv(time = OS.time,event = OS)~',paste0('score',collapse = '+'))%>%as.formula()
train.df$score<-as.numeric(train.df$score)
coxm<-cph(fml,data = train.df,x=T,y=T,surv = T)
med<-Quantile(coxm)
surv<-Survival(coxm)
surv1<-function(x) surv(times = 30*12,lp = x)
surv2<-function(x) surv(times = 30*36,lp = x)
surv3<-function(x) surv(times = 30*60,lp = x)
dd <- datadist(train.df)
options(datadist=dd)
nom <- nomogram(coxm,fun=list(surv1,surv2,surv3),lp=F,funlabel = c("1-year OS","3-year OS",'5-year OS'))
plot(nom,xfrac = 0.3)

####print results
dir.create('output')
write.csv(uni.result,file = 'output/univariate_cox.csv')
