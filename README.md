# EZsinglecell

### Install
```
# Install development version from GitHub:
# install.packages("devtools")
devtools::install_github("chenweng1991/EZsinglecell",auth_token="0da75cd75a04690879e61859bdc38700d5234c55")
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
  - TFfromDBD
- GSEA msigdb
  - msig.db
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
