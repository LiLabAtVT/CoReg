\name{networkFromEdgeList}
\alias{networkFromEdgeList}
\title{Convert an edge list to a graph object}
\description{Function for converting an edge list into a graph object. Input should be a data frame containing two columns. The first column contains the IDs for transcription factors and second the IDs for target genes}
\usage{
networkFromEdgeList(edgeList,simple = T)
}
\arguments{
\item{edgeList}{ Data frame of the network edge list. }
\item{simple}{ a boolean value indicating whether to load the network as simple graph. \code{TRUE} for yes, \code{FALSE} for no. Default value is \code{TRUE}. \code{simple = TRUE} means the self-loops and multiple edges will be removed once the network is loaded. }
}
\details{
This function converts data frame to an igraph graph object. Data frame should contain two columns. The first column contains the IDs for transcription factors and second the IDs for target genes.
}
\value{
An "igraph" object storing the network
}
\seealso{
  \code{\link{CoReg}}
}
\author{Qi Song}
\references{ Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. }
\keyword{convert edge list}
