\name{networkSim}
\alias{networkSim}
\alias{plot.CoReg.rewSim}
\alias{print.CoReg.rewSim}
\alias{summary.CoReg.rewSim}
\title{Rewiring simulation for the network}
\description{This function performs the rewiring simulation for the input}
\usage{
  networkSim(g,nDup,dDup,rewProb,methods=c(),nRep = 1,nThreads =1)
  \method{print}{CoReg.rewSim}()
  \method{summary}{CoReg.rewSim}()
  \method{plot}{CoReg.rewSim}()
}
\arguments{
  \item{g}{ The input transcription network. This should be an "igraph" object returned by the function \code{"loadNetwork"}  }
  \item{nDup}{ Numeric value specifying number of the genes to be duplicated. These duplicated genes will be rewired }
  \item{dDup}{ Numeric value specifying the minimum degree for the genes to be duplicated. This function first ranks the genes based on their degrees. Then randomly selected \code{nDup} genes of which the degrees are greater than dDup. If number of genes above the threshold is smaller than \code{nDup}, all the genes above the threshold will be selected }
  \item{rewProb}{ A vector specifying the a series of rewiring probabilities }
  \item{methods}{ Clustering methods to be compared with CoReg. Available options are: "coregJac" (CoReg + Jaccard index), "coregGeo" (CoReg + geometric index),"coregInv" (CoReg + inverse log-weighted index),"lp" (label propagation),"wt" (walk trap),"eb" (edge betweenness)}
  \item{nRep}{ Number of replicates for run of rewiring simulation }
  \item{nThreads}{ Number of threads for running the simulation. Only valid when "coregGeo" is selected in \code{Methods} }
}
\details{
  This function performs the rewiring simulation for the network. The simulation starts with duplicating \code{nDup} genes with degrees greater than \code{dDup}. Then the duplicated genes will be rewired using a series of rewiring probability specified by \code{rewProb}. The result can be plotted using generic function \code{plot()}
}
\value{
An object of class \code{CoReg.rewSim}
}
\seealso{
  \code{\link{CoReg}}
}
\author{Qi Song}
\references{}
\examples{
data(athNet, package="CoReg")
simRe<-rewSim(athNet,nDup = 50, dDup = 10, c(0.3,0.5),c("coregJac","lp","wt","eb"),2)

# Display simulation result
print(simRe)
simRe$evalResult

# Plot simulation result
plot(simRe)

}
\keyword{network simulation}