% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bubbleplot.R
\name{Generalbubbleplot}
\alias{Generalbubbleplot}
\title{Generalbubbleplot}
\usage{
Generalbubbleplot(pair.data = NULL, cpcol = "celltype2",
  toskip = NULL, genelist = NULL, hlen = 30, llen = 25,
  midlen = 0, limlen = 30, showcmcol = F, titlename = "",
  angl = 0, donormalscale = F, doreturn = F, usefacet = F)
}
\arguments{
\item{pair.data}{the data pair that was prepred by datapair.mk}

\item{cpcol}{The column name for comparison, the default is name}

\item{toskip}{A vector of variables from col(eg,  a cell type), that I don't want to show in bubble plot}

\item{genelist}{A vector of genes that are interesting}

\item{hlen}{The high color span length, default is 30}

\item{llen}{The low color span length, default is 25}

\item{midlen}{The middle color span length, default is 0}

\item{limlen}{The extreme high or low color span length}

\item{showcmcol}{If true, print out a couple of lines of compared column}

\item{titlename}{The figure title}

\item{angl}{The angle of x axis that I want it to rotate}

\item{donormalscale}{If true,  further normalize to zero for all genes}

\item{doreturn}{If true,  reurn a list including bubble data and bubble plot}

\item{usefacet}{If true,  plot with facet, an extra column in genelist is required}
}
\value{
return data pair that can be used for DE, bubble plot and others.
}
\description{
This function is to make bubbleplot out of a data pair.  Notice,  genelist is a dataframe that contain a column that is gene name, where the column name is "gene",  if fa facet plot is wanted , an extra column named "tag" is required
}
\examples{
Generalbubbleplot(ROCKvsnorock.non.paired,cpcol="name",genelist=c("TP53","TNFRSF1A","BAK1","CASP1"),donormalscale=F)
}
