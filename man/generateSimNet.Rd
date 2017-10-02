\name{generateSimNet}
\alias{generateSimNet}
\title{Generate simulated network and test}
\description{ This function generates simulated network with pre-specified modular structure}
\usage{
generateSimNet(mSize,mNum,targetNum,auxNum,prob)
}

\arguments{
  \item{mSize}{Size of module in simulated network(s)}
  \item{mNum}{Number of regulator gene in each module}
  \item{targetNum}{Number of target genes for each regualtor in the module}
  \item{auxNum}{Number of auxiliary genes in simulated network(s)}
  \item{prob}{Co-regulation probability. A numeric value between 0 and 1. Larger value generates network with stronger co-regulation pattern}
}

\details{
  This function first calls generates simulated network with pre-specified modular structure. The properties of network are specified by the arguments. The approach of generating simulated netwoks is a modified version of a previous publication. Modification is made to generate co-regulatory modules in directed network
}

\value{
  A list including following elements:
  \item{el}{Edge list of the simulated network}
  \item{modulePartition}{The pre-specified module labels for regulator genes. They are considered as "ground truth" of co-regulatory mdoules}
}

\seealso{
 \code{\link{netSimAndEval}}
}
\author{Qi Song}
\examples{
re<-generateSimNet(10,5,20,100,0.5)
}