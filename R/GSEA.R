#' Fun.enrich_withFC
#'
#' This function is to do GSEA enrichment analysis
#' @param markergenesList is a list containing several elements each of which is a vector of gene names
#' @param All.genes the background genes  data(all.genes.refer)
#' @param db  =msig.db   the term database to use, default is the msig.db;  data(msig.db)
#' @param qcutoff=0.05  the term to show with q value lower than 0.05
#' @param top=NULL  pass a number to only show top # terms
#' @return this will return a simple object with two elements, one is the table to show the significant terms , the other one shows the gene names involved significant terms
#' @export
#' @examples
#' S7G.GSEA<-Fun.enrich_withFC(markergenesList=list(S7G=topS7G.genes),All.genes=all.genes.refer)
Fun.enrich_withFC<-function(markergenesList,All.genes=all.genes.refer
,db=msig.db,qcutoff=0.05,top=NULL)
{
	require(qvalue)
	binomial.test.resultsList<-list()
	resultname<-c()
	for (markname in names(markergenesList))
	{
		if(length((markergenesList[[markname]]))==0)
		{
			next
		}
		print(markname)
		binomialtest.msig.result<-c()
		for (name in names(db))
		{
			binomialtest.msig.result<-rbind(binomialtest.msig.result,binomialtest.msig.enrch_deplet(as.character(markergenesList[[markname]]),All.genes,name,thedatabase=db))
		}
		fdr.enrich<-qvalue(binomialtest.msig.result[,2])$qvalues
		fdr.deplete<-qvalue(binomialtest.msig.result[,3])$qvalues
		FDR.combined<-c()
		FDR.combined[which(is.na(binomialtest.msig.result$ratio))]<-NA
		FDR.combined[which(binomialtest.msig.result$ratio>1)]<-fdr.enrich[which(binomialtest.msig.result$ratio>1)]
		FDR.combined[which(binomialtest.msig.result$ratio<=1)]<-fdr.deplete[which(binomialtest.msig.result$ratio<=1)]
		if(is.null(FDR.combined))
		{
			next
		}
		binomialtest.msig.result<-cbind(binomialtest.msig.result,FDR=FDR.combined)
		binomialtest.msig.result<-binomialtest.msig.result[order(binomialtest.msig.result$FDR),]
		binomial.test.resultsList<-c(binomial.test.resultsList,list(binomialtest.msig.result))
		resultname<-c(resultname,markname)

	}

	names(binomial.test.resultsList)=resultname
	DRscombo<-data.frame()
	for (markname in names(binomial.test.resultsList))
	{
		tmp<-subset(binomial.test.resultsList[[markname]],FDR<qcutoff)
		row.names(tmp)<-tmp[,"name"]
		colnames(tmp)[c(4,5)]<-paste(markname,colnames(tmp)[c(4,5)],sep="_")
		for (othermarkname in setdiff(names(binomial.test.resultsList),markname))
		{
			tmp2<-binomial.test.resultsList[[othermarkname]]
			row.names(tmp2)<-tmp2[,"name"]
			tmp2<-tmp2[row.names(tmp),c("ratio","FDR")]
			colnames(tmp2)<-paste(othermarkname,colnames(tmp2),sep="_")
			tmp<-cbind(tmp,tmp2)
		}
		tmp<-tmp[,c(-1,-2,-3)]
		if (!is.null(top))
		{
			tmp<-head(tmp,n=top)
		}
		DRscombo<-Tomerge.col(DRscombo,tmp)
	}
	all.list<-c()
	for (term in row.names(DRscombo))
	{
		cur.genelist<-c()
		for (i in 1:length(markergenesList))
		{
			cur.genes<-markergenesList[[i]][markergenesList[[i]] %in% as.character(db[[term]])]
			cur.genelist<-c(cur.genelist,list(cur.genes))
		}
		names(cur.genelist)<-names(markergenesList)
		all.list<-c(all.list,list(cur.genelist))

	}
	names(all.list)<-row.names(DRscombo)
#)
#DRscombo<-DRscombo[tmpcondition,]
	return(list(DRscombo,all.list=all.list))
}



#' Fun.enrich_withFC.pvalue
#'
#' This function is to do GSEA enrichment analysis, to be updated
#' @param markergenesList is a list containing several elements each of which is a vector of gene names
#' @param All.genes the background genes  data(all.genes.refer)
#' @param db  =msig.db   the term database to use, default is the msig.db;  data(msig.db)
#' @param qcutoff=0.05  the term to show with q value lower than 0.05
#' @param top=NULL  pass a number to only show top # terms
#' @return this will return a simple object with two elements, one is the table to show the significant terms , the other one shows the gene names involved significant terms
#' @export
#' @examples
#'
Fun.enrich_withFC.pvalue<-function(markergenesList,All.genes=all.genes.refer,db=msig.db,qcutoff=0.05,top=NULL)
{
	require(qvalue)
	binomial.test.resultsList<-list()
	resultname<-c()
	for (markname in names(markergenesList))
	{
		if(length((markergenesList[[markname]]))==0)
		{
			next
		}
		print(markname)
		binomialtest.msig.result<-c()
		for (name in names(db))
		{
			binomialtest.msig.result<-rbind(binomialtest.msig.result,binomialtest.msig.enrch_deplet(as.character(markergenesList[[markname]]),All.genes,name,thedatabase=db))
		}
		fdr.enrich<-binomialtest.msig.result[,2]
		fdr.deplete<-binomialtest.msig.result[,3]
		FDR.combined<-c()
		FDR.combined[which(is.na(binomialtest.msig.result$ratio))]<-NA
		FDR.combined[which(binomialtest.msig.result$ratio>1)]<-fdr.enrich[which(binomialtest.msig.result$ratio>1)]
		FDR.combined[which(binomialtest.msig.result$ratio<=1)]<-fdr.deplete[which(binomialtest.msig.result$ratio<=1)]
		if(is.null(FDR.combined))
		{
			next
		}
		binomialtest.msig.result<-cbind(binomialtest.msig.result,FDR=FDR.combined)
		binomialtest.msig.result<-binomialtest.msig.result[order(binomialtest.msig.result$FDR),]
		binomial.test.resultsList<-c(binomial.test.resultsList,list(binomialtest.msig.result))
		resultname<-c(resultname,markname)

	}

	names(binomial.test.resultsList)=resultname
	DRscombo<-data.frame()
	for (markname in names(binomial.test.resultsList))
	{
		tmp<-subset(binomial.test.resultsList[[markname]],FDR<qcutoff)
		row.names(tmp)<-tmp[,"name"]
		colnames(tmp)[c(4,5)]<-paste(markname,colnames(tmp)[c(4,5)],sep="_")
		for (othermarkname in setdiff(names(binomial.test.resultsList),markname))
		{
			tmp2<-binomial.test.resultsList[[othermarkname]]
			row.names(tmp2)<-tmp2[,"name"]
			tmp2<-tmp2[row.names(tmp),c("ratio","FDR")]
			colnames(tmp2)<-paste(othermarkname,colnames(tmp2),sep="_")
			tmp<-cbind(tmp,tmp2)
		}
		tmp<-tmp[,c(-1,-2,-3)]
		if (!is.null(top))
		{
			tmp<-head(tmp,n=top)
		}
		DRscombo<-Tomerge.col(DRscombo,tmp)
	}
	all.list<-c()
	for (term in row.names(DRscombo))
	{
		cur.genelist<-c()
		for (i in 1:length(markergenesList))
		{
			cur.genes<-markergenesList[[i]][markergenesList[[i]] %in% as.character(db[[term]])]
			cur.genelist<-c(cur.genelist,list(cur.genes))
		}
		names(cur.genelist)<-names(markergenesList)
		all.list<-c(all.list,list(cur.genelist))

	}
	names(all.list)<-row.names(DRscombo)
#)
#DRscombo<-DRscombo[tmpcondition,]
	return(list(DRscombo,all.list=all.list))
}



#' binomialtest.msig.enrch_deplet
#'
#' This function is an internal function calculating the significance
#' @param mylist
#' @param All=All.genes
#' @param name
#' @param thedatabase=db
#' @return
#' @export
#' @examples
binomialtest.msig.enrch_deplet<-function(mylist,All=All.genes,name,thedatabase=db)
{
n<-length(mylist)
x<-length(which(mylist %in% thedatabase[[name]]))
p<-length(setdiff(thedatabase[[name]],setdiff(thedatabase[[name]],All)))/length(All)

	binomitest.enrich<-binom.test(x,n,p,alternative="greater")
	binomitest.deplete<-binom.test(x,n,p,alternative="less")
statistics<-data.frame(name,enrichsc=binomitest.enrich$p.value,depletesc=binomitest.deplete$p.value,ratio=(x/n)/p)
return(statistics)
}



#' Tomerge.col
#'
#' This function is an internal function doing ther merge
#' @param df1,
#' @param df2
#' @return ""
#' @export
#' @examples
Tomerge.col<-function(df1,df2)
{
	if(length(df1)==0)
	{
		newdf<-df2
		return(newdf)
	}else
	{
		df2<-df2[!row.names(df2) %in%  row.names(df1),]
		newdf<-data.frame(row.names=c(row.names(df1),row.names(df2)))
		for(name in names(df1))
		{
			if(is.factor(df1[,name]))
			{
				assign(name,c(as.character(df1[,name]),as.character(df2[,name])))
			}else
			{
				assign(name,c(df1[,name],df2[,name]))
			}
			newdf<-data.frame(newdf,get(name))
			names(newdf)[ncol(newdf)]<-name
		}
		return(newdf)
	}
}
