##Valiation cohort preprocessing
temp.env <- new.env()
extract.3.genes <- function(){
g15459 <- read.table('validation cohorts/GSE15459_series_matrix.txt',
                     skip = 72,header = T,sep = '\t',encoding = 'UTF-8',fill = T,)%>%.[1:(nrow(.)-1),]
    read.xlsx('validation cohorts/GSE15459_outcome.xls',sheetIndex = 1)
    mx.g15459<-g15459%>%column_to_rownames(var = 'ID_REF')%>%as.matrix()
    # tmm.g15459<-DGEList(mx.g15459)%>%
    #   calcNormFactors(method='TMM')%>%
    #   cpm(log=F)
    genes.cv <- read.table('3genes.txt',header = T,fill = T,sep = '\t')
    gene3.df<-mx.g15459[rownames(mx.g15459)%in%genes.cv$To,]%>%as.data.frame%>%rownames_to_column(var = 'ID_REF')
    gene3.df$ID_REF<-genes.cv$From[match(gene3.df$ID_REF,genes.cv$To)]
    
    genes.df<-aggregate(.~ID_REF,data = gene3.df,FUN = sum)%>%column_to_rownames('ID_REF')%>%as.matrix()%>%t()%>%log1p()
    genes.df%<>%as.data.frame()%>%rownames_to_column(var = 'GSM.ID')
    assign('gene3.df',value = genes.df,envir = temp.env)
}
extract.clin <- function(){
  g15459.clin<-read.xlsx('validation cohorts/GSE15459_outcome.xls',sheetIndex = 1,endRow = 193)
  assign('g15459.clin',value = g15459.clin,envir = temp.env)
}
extract.3.genes()
extract.clin()

assign('g15459.df',merge(temp.env$gene3.df,temp.env$g15459.clin,by = 'GSM.ID'),envir = .GlobalEnv)
rm(temp.env)
