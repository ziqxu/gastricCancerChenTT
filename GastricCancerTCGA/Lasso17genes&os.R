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

wt.coef <- coef(fit,s = cv.fit$lambda.min)%>%as.matrix()%>%as.data.frame%>%filter(`1`!=0)%>%as.matrix()
#update.list <- c("ADAM12","CDH19","SLC5A7","SHOX2","C1QL1","RHPN2","APLN","MNX1","IGF2BP1","HMGB3")
update.list <- rownames(wt.coef)
fml<-paste0('surv~',paste0(update.list,collapse = '+'))%>%as.formula()
#wt.coef<-coef(coxph(fml,data = os.anal.df[train.set,]))%>%as.matrix()
risk.score<-x.mat[,update.list]%*%wt.coef
median(risk.score)

train.df <- os.anal.df[train.set,]
if(all(rownames(train.df)==rownames(risk.score))){
train.df$score <- risk.score[,1]
}
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
  y <- cbind(y$coefficient, y$conf.int,concordance=y$concordance[1],c.se=y$concordance[2],n=y$n)
  y <- y[1:nrow(y),c(2,8,9,5,10,11,12),drop=F]

  return(y)
}

uni.anal('score',train.df)->test
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

mulcox.list <- c('score.group','T_group','N_group','M_group','Age')
mulcox.df<-merge(train.df,pheno, by = 0)%>%
  mutate(M_group = as.factor(pathologic_M)%>%relevel(ref='M0'))%>%
  mutate(T_group = ifelse(grepl('T1',pathologic_T)|grepl('T2',pathologic_T),'≤T2',ifelse(pathologic_T=='TX','TX','>T2')))%>%
  mutate(N_group = ifelse(grepl('N0',pathologic_N),'N0',ifelse(pathologic_N=='NX','NX','>N0')))%>%
  mutate(Age = ifelse(age_at_index.demographic>60,'>60','≤60'))%>%
  mutate(Age=as.factor(Age)%>%relevel(ref='≤60'))%>%
  mutate(T_group=as.factor(T_group)%>%relevel(ref='≤T2'))%>%
  mutate(N_group=as.factor(N_group)%>%relevel(ref='N0'))%>%
  mutate(Sex = as.factor(gender.demographic))%>%
  .[c(mulcox.list,'OS.time','OS')]%>%
  na_if('NX')%>%
  na_if('TX')%>%
  na_if('MX')%>%
  droplevels()
# table(mulcox.df$pathologic_M)
# table(mulcox.df$M_group)
# table(mulcox.df$pathologic_T)
# table(mulcox.df$T_group)
# table(mulcox.df$pathologic_N)
# table(mulcox.df$N_group)


fml<-paste0('Surv(time = OS.time,event = OS)~',paste0(mulcox.list,collapse = '+'))%>%as.formula()
mulcox.df$score<-as.numeric(mulcox.df$score)
coxm<-cph(fml,data = mulcox.df,x=T,y=T,surv = T)
med<-Quantile(coxm)
surv<-Survival(coxm)
surv1<-function(x) surv(times = 30*12,lp = x)
surv2<-function(x) surv(times = 30*36,lp = x)
surv3<-function(x) surv(times = 30*60,lp = x)
dd <- datadist(mulcox.df)
options(datadist=dd)
nom <- nomogram(coxm,fun=list(surv1,surv2,surv3),lp=F,funlabel = c("1-year OS","3-year OS",'5-year OS'))
plot(nom,xfrac = 0.3)

#### DCA ####
require(ggsci)
cols<-ggsci::pal_jama(palette = c("default"))(3)

coxm<-cph(fml,data = mulcox.df,x=T,y=T,surv = T,time.inc=365)
set.seed(seed)
cal<-calibrate(coxm,cmethod = 'KM',method = 'boot',u = 365,m = 100,B = 100)

coxm<-cph(fml,data = mulcox.df,x=T,y=T,surv = T,time.inc =365*3)
set.seed(seed)
cal2<-calibrate(coxm,cmethod = 'KM',method = 'boot',u = 365*3,m = 100,B = 100)

coxm<-cph(fml,data = mulcox.df,x=T,y=T,surv = T,time.inc =365*5)
set.seed(seed)
cal3<-calibrate(coxm,cmethod = 'KM',method = 'boot',u = 365*5,m = 100,B = 100)
pdf(file='output/135_year_DCA.pdf',onefile = T)
{
  plot(cal,lwd=2,lty=1,errbar.col=cols[1],
       xlim=c(0.1,1),ylim=c(0.1,1), legend=F, subtitles=F,
       xlab='Nomogram-predicted OS probability',
       ylab="Actual OS probability",
       col=cols[1]
       )
    lines(cal[,c('mean.predicted','KM')],col=cols[1])
    par(new=TRUE)
    plot(cal2,lwd=2,lty=1,errbar.col=cols[2],
         xlim=c(0.1,1),ylim=c(0.1,1),legend=F, subtitles=F,
         xlab='',
         ylab="",
         col=cols[2]
    )
    lines(cal2[,c('mean.predicted','KM')],col=cols[2])
    par(new=T)
    plot(cal3,lwd=2,lty=1,errbar.col=cols[3],legend=F, subtitles=F,
         xlim=c(0.1,1),ylim=c(0.1,1),
         xlab='',
         ylab="",
         col=cols[3]
    )
    lines(cal3[,c('mean.predicted','KM')],col=cols[3])
    abline(0,1,lty=2,col='gray')
  # lines(cal[,c('mean.predicted','KM')])
    abline(0,1,lty=2)
    legend(0.85,0.3,c(paste("1 year"),
                   paste("3 year"),
                   paste("5 year")),
         x.intersp=1, y.intersp=0.8,
         lty= 1 ,lwd= 2,col=cols,
         bty = "n",# bty框的类型
         seg.len=1,cex=0.8)
}
dev.off()
# pdf(file='output/two_year_DCA.pdf',onefile = T)
# {
#   plot(cal2,lwd=2,lty=1,errbar.col=c(rgb(173,118,192,maxColorValue = 255)),
#        xlim=c(0.1,1),ylim=c(0.1,1),
#        xlab='Nomogram-predition of three-year OS probablity',
#        ylab="Actual OS probability",
#        col=c(192,98,83,maxColorValue=255)
#   )
#   
#   abline(0,1,lty=2)
# }
# dev.off()
# 
# pdf(file='output/5_year_DCA.pdf',onefile = T)
# {
#   plot(cal3,lwd=2,lty=1,errbar.col=c(rgb(56,118,192,maxColorValue = 255)),
#        xlim=c(0.1,1),ylim=c(0.1,1),
#        xlab='Nomogram-predition of 5-year OS probablity',
#        ylab="Actual OS probability",
#        col=c(192,98,83,maxColorValue=255)
#   )
#   lines(cal3[,c('mean.predicted','KM')])
#   abline(0,1,lty=2)
# }
# dev.off()
####cox ####
mul.anal <- function(x,df){
  surv<-Surv(time = df$OS.time,event = df$OS)
  fml<-paste0('surv~',paste0(x,collapse = '+'))%>%as.formula()
  y <- coxph(fml,data = df)%>%summary()
  y <- cbind(y$coefficient, y$conf.int)
  y <- y[,c(1,2,8,9,5),drop=F]
  # y<-coefficients(y)
  return(y)
}
mul.result<-mul.anal(mulcox.list,mulcox.df)
coxph(fml,data = mulcox.df)%>%summary()%>%.$concordance

uni.result<-do.call(rbind,lapply(mulcox.list,function(x)uni.anal(x,mulcox.df)))

#### AUC ####
require(ggsci)
cols<-ggsci::pal_jama(palette = c("default"))(3)

auc.df <- train.df[c('score','OS','OS.time')]
sroc.score<-survivalROC(Stime = train.df$OS.time,status = train.df$OS,marker = train.df$score, predict.time = 365,method = 'KM',)
sroc.score3<-survivalROC(Stime = train.df$OS.time,status = train.df$OS,marker = train.df$score, predict.time = 365*3,method = 'KM',)
sroc.score5<-survivalROC(Stime = train.df$OS.time,status = train.df$OS,marker = train.df$score, predict.time = 365*5,method = 'KM',)

dir.create('output')
pdf('output/TCGA_riskscore_AUC.pdf',onefile = T,width = 6,height = 6)
plot(sroc.score$FP,sroc.score$TP,type="l", xlim=c(0,1), ylim=c(0,1),col=cols[1],   
     xlab=paste( "FP", "\n", "AUC = ",round(sroc.score$AUC,3)), 
     ylab="TP",main="Risk score, Method = KM \n Year = 1")
lines(sroc.score3$FP,sroc.score3$TP,type="l", xlim=c(0,1), ylim=c(0,1),col=cols[2])
lines(sroc.score5$FP,sroc.score5$TP,type="l", xlim=c(0,1), ylim=c(0,1),col=cols[3])
abline(0,1,lty=2,col='gray')
legend(0.6,0.2,c(paste("1 year AUC =",round(sroc.score$AUC,3)),
                 paste("3 year AUC =",round(sroc.score3$AUC,3)),
                 paste("5 year AUC =",round(sroc.score5$AUC,3))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,col=cols,
       bty = "n",# bty框的类型
       seg.len=1,cex=0.8)#

dev.off()
####print results
dir.create('output')
write.csv(uni.result,file = 'output/univariate_cox.csv')
write.csv(mul.result,file = 'output/mulvariate_cox.csv')
pdf(file = 'output/nomogram.pdf',onefile = T,width = 10,height = 7)
plot(nom,xfrac = 0.3)
dev.off()

#### KM-plot
km.plot <- function(df,x){
  {
    fml<-as.formula(paste0('Surv(OS.time,OS)~',x))
    
    data.survdiff <- survdiff(fml, data = df)
    smry = summary(coxph(fml, data = df,method = 'efron'))
    
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    # smry = summary(coxph(fml, data = df,method = 'efron'))
    HR = signif(smry$coef[2], digits=2)
    up95 = smry$conf.int[,"upper .95"]
    low95 = smry$conf.int[,"lower .95"]
    HR <- paste("HR:", round(HR,2), sep = "")
    CI <- paste("95%CI: ", paste(round(low95,2),round(up95,2), sep = ifelse(low95>1,'*-','-')), sep = "")
    HRCI <- paste0(HR,"(",CI,")")
    colorset<-c('#0909a0','#7a0505')
    
  }
  p <- ggsurvplot(fit = surv_fit(fml,data = df),
                  palette = colorset,
                  pval=T,
                  #pval=paste(HRCI,pval = ifelse(p.val < 0.001, "p < 0.001", paste("p = ",round(p.val,4), sep = "")), sep = "\n"),
                  pval.method = T,
                  #    pval.coord=c(max(na.omit(df$OS_MONTHS))*0.6,0.8),
                  risk.table = T,
                  title=paste0('KM-plot of Risk Score')
  )+labs(x='Days')
  
  p$table <- p$table + theme(axis.line = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank(),
                             axis.text.x = element_blank(),
  )
  
  p
}
pdf(file='output/kmplot_tcga.pdf',onefile = T,width = 8,height = 6)
print(km.plot(mulcox.df,x = 'score.group'))
dev.off()

#### 
save(wt.coef,file = 'weight.coeff.RData')
