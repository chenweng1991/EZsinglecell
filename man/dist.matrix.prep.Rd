% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TreCCAtree.R
\name{dist.matrix.prep}
\alias{dist.matrix.prep}
\title{dist.matrix.prep}
\usage{
dist.matrix.prep(binary.primary, datainfo.col1 = c("res.0.6", "nUMI",
  "Sample"), cluster.col1 = "res.0.6", datainfo.col2 = c("res.0.6",
  "nUMI", "Sample"), cluster.col2 = "res.0.6", res1 = "res.0.06",
  res2 = "res.0.06", datainfo.col3 = c("res.0.6", "nUMI", "Sample"))
}
\arguments{
\item{binary.primary}{This is the primary clustering result that can be retrieved from Tree.build.prepare object. (Because Tree.build.prepare mainly do the primary clustering.)The input here is tree.prep.ob$primary.total}

\item{datainfo.col1}{default is c("res.0.6","nUMI","Sample")}

\item{cluster.col1}{default is "res.0.6"}

\item{datainfo.col2}{default is c("res.0.6","nUMI","Sample")}

\item{cluster.col2}{default is "res.0.6"}

\item{res1}{default is "res.0.06}

\item{res2}{default is "res.0.06"}

\item{datainfo.col3}{c("res.0.6","nUMI","Sample")}
}
\value{
This will return an object used for Tree.build.1
}
\description{
This is an important internal function to calculate distance matrix, it is used in Tree.build.1 and Tree.build.2nd.clustering
}
\examples{
 	total.dist.matrxes<-dist.matrix.prep(primary.total,datainfo.col1=c(paste("res.",first.reso[1],sep=""),"nUMI","Sample"),cluster.col1=paste("res.",first.reso[1],sep=""),datainfo.col2=c(paste("res.",first.reso[2],sep=""),"nUMI","Sample"),cluster.col2=paste("res.",first.reso[2],sep=""),res1=paste("res.",first.reso[1],sep=""),res2=paste("res.",first.reso[2],sep=""))
}
