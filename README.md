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

1. To be continue


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
       
### All functions included 
- Fullplot_v2             
- GetPCAcelldata_v2       
- GetPCheat              
- Gettopgenes             
- GettsnesignatureSuper   
- Tomerge                 
- Tomerge_v2              
- format2 


