### Pre-processing of TCGA htseq count data

while(!require(edgeR)){install.packages('edgeR')}
while(!require(magrittr)){install.packages('magrittr')}
while(!require(biomaRt)){install.packages('biomaRt')}
while(!require(tibble)){install.packages('tibble')}
while(!require(dplyr)){install.packages('dplyr')}

raw.count <- read.table('TCGA-STAD.htseq_counts.tsv',header = T)
raw.count.mx <- raw.count%>%column_to_rownames(var="Ensembl_ID")%>%as.matrix()
raw.count.mx.reverted <- (2^raw.count.mx)-1
tmm.mx<-DGEList(raw.count.mx.reverted)%>%
  calcNormFactors(method='TMM')%>%
  cpm(log=T)

ensg.id<-rownames(tmm.mx)%>%sub(pattern = '\\..*','',.)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = ensg.id, mart = mart)

gene_IDs$hgnc_symbol<-na_if(gene_IDs$hgnc_symbol,y = '')
gene_IDs$hgnc_symbol[duplicated(gene_IDs$hgnc_symbol)&&!is.na(gene_IDs$hgnc_symbol)]
no.dup.ensembl.id<-!duplicated(gene_IDs$ensembl_gene_id)
gene_IDs <- gene_IDs[no.dup.ensembl.id,]

tmm.df<-tmm.mx%>%as.data.frame()
tmm.df$ensg.id<- ensg.id
tmm.hugo.df <- merge(gene_IDs,tmm.df,by.y = 'ensg.id',by.x='ensembl_gene_id')
table(duplicated(tmm.hugo.df$ensembl_gene_id))

save(tmm.hugo.df,file = "gastric_tcga_tmm.RData")
