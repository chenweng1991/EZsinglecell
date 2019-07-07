## Loading packages
```
library(EZsinglecell)
library(dplyr)
library(biomaRt)
library(gage)
```
## Prepare insulin regulator gene by Crispr-screening
```
Crisp.t1<-read.table("BInomial.tier1")
Crisp.t2<-read.table("BInomial.tier2")
crisperlist.dic<-read.table("data.cnsis.report")
crisperlist.dic<-crisperlist.dic[,c(1,5)]
crisperlist.dic<-crisperlist.dic[!duplicated(crisperlist.dic$gene_id),]
crisperlist.dic<-data.frame(crisperlist.dic$Direction,row.names=crisperlist.dic$gene_id)
Crisp.t1<-Tomerge_v2(Crisp.t1,crisperlist.dic)
Crisp.t2<-Tomerge_v2(Crisp.t2,crisperlist.dic)
```
## Prepare DB/OB related GWAS data
```
dups<-function(df,col=1){
df<-df[!duplicated(df[,col]),]
row.names(df)<-df[,col]
df<-df[,-col]
return(df)
}
#---
GWASdata<-dups(read.table("GWAS2717.geneinformation",sep="\t",fill=T,header=F))
colnames(GWASdata)<-c("GWAS_cat","GWAS_cat2")
```
## Prepare surfaceome data
```
Surfaceome.data<-as.data.frame(read.table("/mnt/NFS/homeGene/JinLab/cxw486/DBs/Surfasome/Surfaceome.data.csv")[,1])
row.names(Surfaceome.data)<-Surfaceome.data[,1]
names(Surfaceome.data)<-"surface.genes"
```
## Prepare cell cycle data
```
cellcycles<-read.csv("cellcycle.csv",header=T)
G1.S<-gsub(" ","", as.character(cellcycles[,1])) %>% .[nchar(.)>0]
S<-gsub(" ","", as.character(cellcycles[,2])) %>% .[nchar(.)>0]
G2.M<-gsub(" ","", as.character(cellcycles[,3])) %>% .[nchar(.)>0]
M<-gsub(" ","", as.character(cellcycles[,4])) %>% .[nchar(.)>0]
M.G1<-gsub(" ","", as.character(cellcycles[,5])) %>% .[nchar(.)>0]
all.cellcycle<-unique(c(G1.S,S,G2.M,M,M.G1))
```
## Prepare Transcription factor
```
TFdatabase<-read.table("~/DBs/DBD_transcriptionfactorDB/update.18.12.18/hs.tf.ass")
TFdatabase.2<-read.delim("~/DBs/DBD_transcriptionfactorDB/update.18.12.18/Homo_sapiens_TF.txt")
TFdatabase.3<-getBM(attributes= "hgnc_symbol",
filters=c("go"),
values="GO:0003700", mart=mart)
TFfromDBD<-TFfromDBD[!duplicated(TFfromDBD$hgnc_symbol),]
row.names(TFfromDBD)<-TFfromDBD$hgnc_symbol
TFvector<-unique(c(as.character(TFfromDBD$hgnc_symbol),levels(TFdatabase.2$Symbol),TFdatabase.3$hgnc_symbol))
```
## Prepare GSEA msigdb
```
msig.db<-readList("~/DBs/MsigDB/msigdb.v5.2.symbols.gmt")
```
## Prepare Beta/alpha/delta/pp/episilon specific geneset
```
all.Beta.maker<-readRDS("all.Beta.maker.RDS")
allAlpha.marker<-readRDS("allAlpha.marker")
delta.markers<-readRDS("delta.markers")
pp.markers<-readRDS("pp.markers")
epsilon.markers<-readRDS("epsilon.markers")
```
## Prepare Diabetes and obesity trajectory genes
```
Beta.BMI.trj.genes<-read.csv("Beta.BMI.trj.genes.csv",row.names=1)
Beta.T2D.trj.genes<-read.csv("Beta.T2D.trj.genes.csv",row.names=1)
```
## Prepare druggable gene
```
gene_drug.db<-read.delim("interactions.tsv",header=T)
```
## Prepare hocomoco genes
```
Hocomocogenes<-read.delim("HUMAN_mono_motifs.csv")
```
## prepare GSEA all gene refer
ref.gene.tb<-read.delim("/mnt/NFS/homeGene/JinLab/cxw486/Genome/hg18/refgene/briefGene.table")
all.genes.refer<-unique(ref.gene.tb$name2)



# Save data
devtools::use_data(Crisp.t1,Crisp.t2,GWASdata,Surfaceome.data,G1.S,S,G2.M,M,M.G1,all.cellcycle,TFfromDBD,msig.db,all.Beta.maker,allAlpha.marker,delta.markers,pp.markers,epsilon.markers,Beta.BMI.trj.genes,Beta.T2D.trj.genes,gene_drug.db,Hocomocogenes,overwrite=T)
devtools::use_data(TFvector)
devtools::use_data(all.genes.refer)
