source('load_package.R')
load('weight.coeff.RData')
gene.list <- rownames(wt.coef)
####load gse14359####
source('extract_g15459.R')



####cox.df####
cox.df <- g15459.df%>%
  select(all_of(gene.list), OS = Outcome..1.dead., OS.time= Overall.Survival..Months...)
{
  rep.df<-cox.df
  groupby_median <- function(x){
    median_x <- median(x)
    result <- ifelse(x<median_x,0,1)
    return(result)
  }
  rep.df[gene.list] <- lapply(rep.df[gene.list],groupby_median)
  cox.df<-rep.df
  rm(rep.df)
  }
cox.df$risk.score <- as.matrix(cox.df[gene.list])%*%wt.coef%>%.[,1]
cox.df$group <- ifelse(cox.df$risk.score>median(cox.df$risk.score),1,0)
#####COX#####
uni.anal <- function(x,df){
  fmt.str <- paste0('Surv(time = OS.time,event = OS)~',x)%>%
    as.formula()
  y <- coxph(fmt.str, data = df)%>%summary()
  y <- cbind(y$coefficient, y$conf.int,concordance=y$concordance[1],c.se=y$concordance[2],n=y$n)
  y <- y[1:nrow(y),c(2,8,9,5,10,11,12),drop=F]
  
  return(y)
}
uni.anal('group',df = cox.df)

#####KMplot######
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
km.plot(df = cox.df,x = 'group')
km.plot(df = cox.df,'SHOX2')

####ROC AUC####
require(ggsci)
cols<-ggsci::pal_jama(palette = c("default"))(3)
dd <- datadist(cox.df)
options(datadist=dd)

auc.df <- cox.df[c('risk.score','OS','OS.time')]
sroc.score<-survivalROC(Stime = cox.df$OS.time,status = cox.df$OS,marker = cox.df$risk.score, predict.time = 12,method = 'KM',)
sroc.score3<-survivalROC(Stime = cox.df$OS.time,status = cox.df$OS,marker = cox.df$risk.score, predict.time = 12*3,method = 'KM',)
sroc.score5<-survivalROC(Stime = cox.df$OS.time,status = cox.df$OS,marker = cox.df$risk.score, predict.time = 12*5,method = 'KM',)

dir.create('output')
#pdf('output/TCGA_riskscore_AUC.pdf',onefile = T,width = 6,height = 6)
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

#dev.off()

