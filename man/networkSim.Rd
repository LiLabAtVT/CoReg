\name{networkSim}
\alias{networkSim}
\alias{plot.CoReg.rewSim}
\alias{print.CoReg.rewSim}
\alias{getRewScore.CoReg.rewSim}
\title{Rewiring simulation for the network}
\description{This function performs the rewiring simulation for the input}
\usage{
  networkSim(g,nDup,dDup,rewProb,compareMethods=c(),nRep = 1,nThreads =1)
}
\arguments{
  \item{g}{ The input transcription network. This should be an "igraph" object returned by the function \code{"loadNetwork"}  }
  \item{nDup}{ numeric value specifying number of the genes to be duplicated. These duplicated genes will be rewired }
  \item{dDup}{ numeric value specifying the minimum degree for the genes to be duplicated. This function first ranks the genes based on their degrees. Then randomly selected \code{nDup} genes of which the degrees are greater than dDup. If number of genes above the threshold is smaller than \code{nDup}, all the genes above the threshold will be selected }
  \item{rewProb}{ a vector specifying the a series of rewiring probabilities }
  \item{compareMethods}{ clustering methods to be compared with CoReg. Available options are: "lp" (label propagation),"wt" (walk trap),"eb" (edge betweenness)}
  \item{nRep}{ number of replicates for run of rewiring simulation }
  \item{nThreads}{ number of threads for running the simulation }
}
\details{
  This function performs the rewiring simulation for the network. The simulation starts with duplicating \code{nDup} genes with degrees greater than \code{dDup}. Then the duplicated genes will be rewired using a series of rewiring probability specified by \code{rewProb}. The result can be plotted using generic function \code{plot()}
}
\value{
An object of class \code{CoReg.rewSim}, a list including elements
\item{plotData}{ a dataframe storing the data for plotting the result }
\item{evalResult}{ a dataframe reporting the simulation result }
}
\seealso{
  \code{\link{CoReg}}
}
\author{Qi Song}
\references{}
\examples{
data(athNet, package="CoReg")
simRe<-networkSim(athNet,nDup = 50, dDup = 10, c(0.3,0.5),c("lp","wt","eb"),2)

# Display simulation result
print(simRe)
simRe$evalResult

# Plot simulation result
plot(simRe)

}
\keyword{network simulation}