source('load_package.R')
load('weight.coeff.RData')
##Valiation 
g13861<-read.table('validation cohorts/GSE13861_series_matrix.txt',
           skip = 67,header = T,sep = '\t',encoding = 'UTF-8',fill = T,)%>%.[1:(nrow(.)-1),]
g13861.clin <- read.xlsx('validation cohorts/GSE13861_GE_MDACC_DepSysB_ClinicalInformation_update_ver2.xls',
                         sheetIndex = 2,header = T,startRow = 1)
g15459.clin<-read.xlsx('validation cohorts/GSE15459_outcome.xls',sheetIndex = 1)
g15459 <- read.table('validation cohorts/GSE15459_series_matrix.txt',
                     skip = 72,header = T,sep = '\t',encoding = 'UTF-8',fill = T,)%>%.[1:(nrow(.)-1),]
#write.table(g15459$ID_REF,file = 'affy_genelist.txt',sep = '\t',quote = T)

genes.cv <- read.table('3genes.txt',header = T,fill = T,sep = '\t') ##The files need to custom by my self, credit to DAVID gene converter
#g15459%<>%filter(ID_REF%in%genes.cv$To)

g15459$ID_REF<-genes.cv$From[match(g15459$ID_REF,genes.cv$To)]
gene3.df$ID_REF<-genes.cv$From[match(gene3.df$ID_REF,genes.cv$To)]
#dplyr::group_by(g15459,ID_REF)%>%sum()
#dplyr::group_by(gene3.df,ID_REF)%>%sum()

test<-aggregate(.~ID_REF,data = g15459,FUN = sum)%>%column_to_rownames('ID_REF')%>%as.matrix()%>%t()
test<-aggregate(.~ID_REF,data = gene3.df,FUN = sum)%>%column_to_rownames('ID_REF')%>%as.matrix()%>%t()%>%log1p()

read.xlsx('validation cohorts/GSE15459_outcome.xls',sheetIndex = 1)





