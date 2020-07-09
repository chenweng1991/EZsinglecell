# EZsinglecell



### Install
```
# Install development version from GitHub:
# install.packages("devtools")
devtools::install_github("chenweng1991/EZsinglecell")
other requirement:  raster

```
---
### For general convenience
1. Merge dataframe by row.names, only mached rows are shown
```
Tomerge(A,B)
```
2. Merge dataframe by row.names, show all dataframe A rows.
```
Tomerge_v2(A,B)
```
3. Modify plot from ggplot into a clean and squared version
```
format2(p, data, x, y, nolegend = T)
```

### Clustering

1.  Do clustering for one dge sample
```
docluster(dgepreprocess(s7.RockII_1.dge,500,norowname=T),
GetinformativeGene(dgepreprocess(s7.RockII_1.dge,500,norowname=T),500),
"s7.RockII_1",reso=0.6)
```

2.  Do clustering for combined multiple Samples, Number is informative gene number.
```
docluster.multi(Number=500,txcutoff=500,sets=list(s7.RockII_1=s7.RockII_1.dge,s7.B=s7.B.dge),nms=c("s7.RockII_1","s7.B"),israw=T)
```

### DE analysis
1. Make data pair that can be used for DE analysis and bubble plot
```
ROCKvsnorock.endo.paired<-datapair.mk(list(S7rock=S7rock_1.ob,S7=S7.ob),cols=c("Sample","Sample.2nd"),pick.list=list(c("s7.RockII_1"),c("s7.B_1")),normalizecellsize=F,randomizecelloirder=T)
```
2. Run negtive binomial based differentially expression gene calling
```
library(scran)  # Use computeSumFactors to compute size factor
ROCKvsnorock.endo.tri.dummy<-DE.gettripple(ROCKvsnorock.endo.paired,cpcol="name")
ROCKvsnorock.endo.de<-DoDE(ROCKvsnorock.endo.tri.dummy,"name",onlyoneSample=T,cpus=16)
```

### Ploting

1. Fullplot_v2,  plot basic clustering results for seurat object
```
mytopgenes<-Gettopgenes(object, number)  # Sometimes, precalculated cluster gene signatures can save some time
Fullplot_v2(S7rock_1.ob,"S7rock_1.ob.pdf",topgene=NULL,resolusion="res.0.6",signiture=c("INS","GCG","SST","PPY","KRT19","COL1A2","REG1A","DNAJB1","GHRL"))
```
2.  GettsnesignatureSuper, plot tsne or PCA stained by gene relative expression
```
GettsnesignatureSuper(object, object.all, signiture = c("INS", "RBP4","FFAR4", "ID1", "ID2", "ID3", "DNAJB1"), doPCA = F, dotSNE = T,usePC34 = F, extratitle = "", toreturn = F, dotsize = 0.3,buttomgrey = T, nolegend = T, highcolor = "red")
```

3. Generalbubbleplot,  to plot bubble plot for data pair, df is a dataframe whose row name is gene name, there is an extra columns named "tag", that is used for facet if usefacet is set as T.
```
Generalbubbleplot(ROCKvsnorock.non.paired,cpcol="name",genelist=df,donormalscale=F,usefacet=F)
```

### RePACT

0. Loading packages
```
library(pscl)  # For Repact analysis to calculate logistic regression pseudo-p value
library(plot3D) # To plot 3D
```
1. Prepareforpseudoregress.g  This function is to perform the initial regression to prepare the trajectory study
```
T2D.tjct.ob<-Prepareforpseudoregress.g(Beta.HSnegonly.ob,PCrange=1:10,phenodic.use=phenodic,pheno="disease",linear=F)
```

2. Toplot3Dtjct   his function is to make 3D example plot for regression analysis
```
Toplot3Dtjct(T2D.tjct.ob,PCrange=c(1,3,4),pheno="disease",linear=F,multiplotname="test.pdf",titlename="tittle")
```


3. Tjct.core.gen  This function is to compute significant trajectory genes by linear regression
```
T2D.tjct.2nd.ob<-Tjct.core.gen(T2D.tjct.ob)
```

4. Tjct.core.plot  This function is to generate major plots for Repact analysis
```
Tjct.core.plot(T2D.tjct.ob,T2D.tjct.2nd.ob,pheno="Disease",f1.name="T2D.tjct.10d.violin.pdf",f2.name="T2D.tjct.his.pdf",f3.name="T2D.tjct.trj.heatmap.pdf",f3.height=14,f3.tittle="cell type:Changing genes on phenotype trajectory\ntop6%",table1.name="T2D.tjct.traj.up.genes-q0.05Full.csv",table2.name="T2D.tjct.traj.dowb.genes-q0.05Full.csv",rankcut=0.04)
```


### GSEA analysis (visualization functions not included yet)

1. The input is a list of several elements each of which contain a many gene names
```
data(all.genes.refer)
data(msig.db)
S7G.GSEA<-Fun.enrich_withFC(markergenesList=list(S7G=topS7G.genes),All.genes=all.genes.refer)
```

### TreCCA tree
1. Analyze cluster-cluster relationships between any adjacent time point (Use S6.AD and S7.B as example)
  - **1-1 Low resolution clustering and connection analysis(To obtain two objects: *.tree.prep and *.tree.1.ob)**   *Input: dge1, dge2, appropriate resolution*
```
S6.AD_S7.B.tree.prep<-Tree.build.prepare(dge1=S6.AD.dge,dge2=s7.B.dge,name1="s6.AD",name2="s7.B",first.reso=c(0.03,0.03))
S6.AD_S7.B.tree.1.ob<-Tree.build.1(S6.AD_S7.B.tree.prep)
```
Above two line will generate a folder and within the folder, thare are 3 PDF files:Clustering results for each time point and 1 connection result.

  - **1-2 High resolution clustering and connection analysis(To obtain two objects: *.list.patched and *.ob)**   *Input:.tree.prep, .tree.1.ob from Step1*
```
S6.AD_S7.B.tree.2nd.primary_0.3_0.3.list<-Tree.build.2nd.clustering(S6.AD_S7.B.tree.prep,S6.AD_S7.B.tree.1.ob,second.reso=c(0.3,0.3))
S6.AD_S7.B.tree.2nd.primary_0.3_0.3.list.patched<-Tree.build.2nd.clustering.patch(S6.AD_S7.B.tree.2nd.primary_0.3_0.3.list,S6.AD_S7.B.tree.prep,S6.AD_S7.B.tree.1.ob,singledog.reso=0.06)
S6.AD_S7.B.tree.2nd.treemade_0.3_0.3.ob<-Tree.build.2nd.treemaking(S6.AD_S7.B.tree.2nd.primary_0.3_0.3.list.patched,second.reso=c(0.3,0.3),upstremename="s6.AD",downstremename="s7.B",dir="")
```
Note: 1. upstremename=name1 ,downstremename=s7.B; 2. Dir: keep consistent with Tree.build.prepare, where the Child folder was generated at Step 1

2. Create a list object containing dges for each clusters across all stages
  - Load all high resolution *.list.patched
```
H1_S0.tree.2nd.primary_0.06.list.patched<-readRDS("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/Entrance/Esderived/Esdrived_DGE/AttemptFrom17.8.28/RDS/H1_S0.tree.2nd.primary_0.06.list.patched.rds")
S0_S1.24.tree.2nd.primary_0.06.list.patched<-readRDS("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/Entrance/Esderived/Esdrived_DGE/AttemptFrom17.8.28/RDS/S0_S1.24.tree.2nd.primary_0.06.list.patched.rds")
S1.24_S1.tree.2nd.primary_0.06.list.patched<-readRDS("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/Entrance/Esderived/Esdrived_DGE/AttemptFrom17.8.28/RDS/S1.24_S1.tree.2nd.primary_0.06.list.patched.rds")  # S1_D2 with 1 cluster
S1_S2.24.tree.2nd.primary_0.06.list.patched<-readRDS("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/Entrance/Esderived/Esdrived_DGE/AttemptFrom17.8.28/RDS/S1_S2.24.tree.2nd.primary_0.06.list.patched.rds")  # S1_D2 with 1 cluster
S2.24_S2.48.tree.2nd.primary_0.06_0.03.list.patched<-readRDS("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/Entrance/Esderived/Esdrived_DGE/AttemptFrom17.8.28/RDS/S2.24_S2.48.tree.2nd.primary_0.06_0.03.list.patched.RDS")
S2.48_S2.all.tree.2nd.primary_0.06.list.patched<-readRDS("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/Entrance/Esderived/Esdrived_DGE/AttemptFrom17.8.28/RDS/S2.48_S2.all.tree.2nd.primary_0.06.list.patched.rds")
S2.all_S3.all.tree.2nd.primary_0.06.list.patched<-readRDS("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/Entrance/Esderived/Esdrived_DGE/AttemptFrom17.8.28/RDS/S2.all_S3.all.tree.2nd.primary_0.06.list.patched.rds")
S3.all_S4.B.tree.2nd.primary_0.06_0.3.list.patched<-readRDS("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/Entrance/Esderived/Esdrived_DGE/AttemptFrom17.8.28/RDS/S3.all_S4.B.tree.2nd.primary_0.06_0.3.list.patched.RDS")
#I made a mistake which makes S4.B_S5.all.tree.2nd.primary_0.06_0.3.list.patched not accessible S4.B_S5.all.tree.2nd.primary_0.06_0.3.list.patched<-readRDS("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/Entrance/Esderived/Esdrived_DGE/AttemptFrom17.8.28/RDS/S4.B_S5.all.tree.2nd.primary_0.06_0.3.list.patched.RDS")
S5.all_S6.A.tree.2nd.primary_0.3_0.3.list.patched<-readRDS("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/Entrance/Esderived/Esdrived_DGE/AttemptFrom17.8.28/RDS/S5.all_S6.A.tree.2nd.primary_0.3_0.3.list.patched.RDS")
S6.A_S7.B.tree.2nd.primary_0.3_0.3.list.patched<-readRDS("/mnt/NFS/homeGene/JinLab/cxw486/Dropseq/Entrance/Esderived/Esdrived_DGE/AttemptFrom17.8.28/RDS/S6.A_S7.B.tree.2nd.primary_0.3_0.3.list.patched")
```

```
rawdata.lst.hi<-list(H1_S0.tree.2nd.primary_0.06.list.patched$H1_0_s0_0$Seurat.list$H1_0@raw.data)
totalobjects.lst.hi<-list(H1_S0.tree.2nd.primary_0.06.list.patched$H1_0_s0_0$Seurat.list$H1_0)
clusternameinfo.hi<-H1_S0.tree.2nd.primary_0.06.list.patched$H1_0_s0_0$Seurat.list$H1_0@data.info %>% cbind(.,Sample.2nd=paste(.$Sample,.[,5],sep="_")) %>% .$Sample.2nd %>% as.character() %>% unique
Samplenameinfo.hi<-H1_S0.tree.2nd.primary_0.06.list.patched$H1_0_s0_0$Seurat.list$H1_0@data.info$Sample %>% unique
n=1
for (prep in list(H1_S0.tree.2nd.primary_0.06.list.patched,S0_S1.24.tree.2nd.primary_0.06.list.patched,S1.24_S1.tree.2nd.primary_0.08.list.patched,S1_S2.24.tree.2nd.primary_0.08.list.patched,S2.24_S2.48.tree.2nd.primary_0.06_0.03.list.patched,S2.48_S2.all.tree.2nd.primary_0.06.list.patched,S2.all_S3.all.tree.2nd.primary_0.06.list.patched,S3.all_S4.B.tree.2nd.primary_0.06_0.3.list.patched,S5.all_S6.A.tree.2nd.primary_0.3_0.3.list.patched,S6.A_S7.B.tree.2nd.primary_0.3_0.3.list.patched))
{
print(n)
n=n+1
for (thread in names(prep)[grepl("singledog.singledog.seurat.ob",names(prep))])
	{
		totalobjects.lst.hi<-c(totalobjects.lst.hi,list(prep[[thread]]))
		cur.datainfo<-cbind(prep[[thread]]@data.info,Sample.2nd=paste(prep[[thread]]@data.info$Sample,prep[[thread]]@data.info[,5],sep="_"))
		Samplenameinfo.hi<-c(Samplenameinfo.hi,cur.datainfo$Sample %>% unique)
		for (cluster in levels(cur.datainfo$Sample.2nd))
		{
			rawdata.lst.hi<-c(rawdata.lst.hi,list(subset(cur.datainfo,Sample.2nd==cluster) %>% row.names %>% prep[[thread]]@raw.data[,.]))
			clusternameinfo.hi<-c(clusternameinfo.hi,cluster)
		}
	}
for (thread in names(prep)[!grepl("singledog",names(prep))])
	{
		totalobjects.lst.hi<-c(totalobjects.lst.hi,list(prep[[thread]]$Seurat.list[[2]]))
		cur.datainfo<-cbind(prep[[thread]]$Seurat.list[[2]]@data.info,Sample.2nd=paste(prep[[thread]]$Seurat.list[[2]]@data.info$Sample,prep[[thread]]$Seurat.list[[2]]@data.info[,5],sep="_"))
		Samplenameinfo.hi<-c(Samplenameinfo.hi,cur.datainfo$Sample %>% unique)
		for (cluster in levels(cur.datainfo$Sample.2nd))
		{
			rawdata.lst.hi<-c(rawdata.lst.hi,list(subset(cur.datainfo,Sample.2nd==cluster) %>% row.names %>% prep[[thread]]$Seurat.list[[2]]@raw.data[,.]))
			clusternameinfo.hi<-c(clusternameinfo.hi,cluster)
		}
	}
}
```


### Data included
- Insulin regulator gene by Crispr-screening
  - Crisp.t1
  - Crisp.t2
- DB/OB related GWAS data
  - GWASdata
- Surfaceome data
  - Surfaceome.data
- Prepare cell cycle data
  - G1.S
  - S
  - G2.M
  - M
  - M.G1
  - all.cellcycle
- Transcription factor
  - TFfromDBD(Old version)
  - TFvector(version 18.12.18)
- GSEA msigdb
  - msig.db
  - all.genes.refer
- Beta/alpha/delta/pp/episilon specific geneset
  - all.Beta.maker
  - allAlpha.marker
  - delta.markers
  - pp.markers
  - epsilon.markers
- Diabetes and obesity trajectory genes
  - Beta.BMI.trj.genes
  - Beta.T2D.trj.genes
- Druggable gene
  - gene_drug.db
- Hocomoco genes
  - Hocomocogenes
### All functions included
- Fullplot_v2
- GetPCAcelldata_v2
- GetPCheat
- GetinformativeGene
- Gettopgenes
- GettsnesignatureSuper
- Tomerge
- Tomerge_v2
- dgepreprocess
- docluster
- docluster.multi
- format2
- datapair.mk
- DE.gettripple
- DoDE

### Other home made functions

- RunHomeGSEA (Input a ranked dataframe, and another vector for enrichment analysis, out put athe GSEA-like plot and permutatiohn-based p-value)
```
#RankedGeneDF is a dataframe with a column"gene", and ranked by another column "score", 
#TestGenes  is a vector of genes that is a subset of the "gene" column of the RankedGeneDF

RunHomeGSEA<-function(RankedGeneDF,TestGenes,PermutationN=1000){
Step.Positive<-sqrt((nrow(RankedGeneDF)-length(TestGenes))/length(TestGenes))
Step.Negative<-sqrt(length(TestGenes)/(nrow(RankedGeneDF)-length(TestGenes)))
ES.cumulative<-c()
StepAction<-c()
StepAction[which(RankedGeneDF$gene %in% TestGenes)]<-Step.Positive
StepAction[which(!RankedGeneDF$gene %in% TestGenes)]<--Step.Negative
ES.cumulative<-cumsum(c(0,StepAction[-length(StepAction)]))
plot<-data.frame(x=1:length(ES.cumulative),ES=ES.cumulative) %>% ggplot(.)+aes(x,ES)+geom_line()+ggtitle(cluster)
ES.max<-max(ES.cumulative)
## Compute nulldistributation
random.ES.max<-c()
for(i in 1:PermutationN){
randomAllGene<-as.character(sample(RankedGeneDF$gene))
RD.ES.cumulative<-c()
StepAction<-c()
StepAction[which(randomAllGene %in% TestGenes)]<-Step.Positive
StepAction[which(!randomAllGene %in% TestGenes)]<--Step.Negative
RD.ES.cumulative<-cumsum(c(0,StepAction[-length(StepAction)]))
random.ES.max<-c(random.ES.max,max(RD.ES.cumulative))
}
cur.pvalue<-1-length(which(random.ES.max<ES.max))/length(random.ES.max)
return(list(plot=plot,ES=ES.max,pvalue=cur.pvalue))
}
```

