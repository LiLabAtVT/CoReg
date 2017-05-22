\name{computeAuROC}
\alias{computeAuROC}
\alias{print.CoReg.auROC}
\alias{plot.CoReg.auROC}
\alias{getAuc.CoReg.auROC}
\title{compute AUC and ROC for the network}
\description{ This function performs the same rewiring simulation as \code{networkSim}. But computes the AUC and ROC as the result}
\usage{
computeAuROC(g,nDup,dDup,rewProb,simMethods,nThreads=1)
  \method{print}{CoReg.auROC}()
  \method{plot}{CoReg.auROC}()
  \method{getAuc}{CoReg.auROC}()
}
\arguments{
  \item{g}{ The input transcription network. This should be an "igraph" object returned by the function \code{"loadNetwork"}  }
  \item{nDup}{ numeric value specifying number of the genes to be duplicated. These duplicated genes will be rewired }
  \item{dDup}{ numeric value specifying the minimum degree for the genes to be duplicated. This function first ranks the genes based on their degrees. Then randomly selected \code{nDup} genes of which the degrees are greater than dDup. If number of genes above the threshold is smaller than \code{nDup}, all the genes above the threshold will be selected }
  \item{rewProb}{ a numeric value specifying the rewiring probability }
  \item{simMethods}{ similarity indices used. Available options are: "jaccard" (jaccard similarity index),"geometric" (geometric similarity index),"invlogweighted" (inverse log weighted similarity index), "wt" (walk-trap-based similarity index)}
  \item{nThreads}{ number of threads for running the simulation. Only valid when \code{simMethods = "geometric"} }
}
\details{
  This function performes the rewring simulation similar to \code{networkSim}, but uses the AUC values and ROC curves as the evaluation metrics. ROC curves can be plotted using the generic function \code{plot()}. 
}
\value{
An object of class \code{CoReg.auROC}, a list including elements
\item{curve}{ the computed curve for ROC lot }
\item{AUC}{ AUC values }
}
\seealso{
  \code{\link{networkSim}}
}
\author{Qi Song}
\references{ }
\examples{
data(athNet)
re <- computeAuROC(athNet,nDup=50,dDup=10,rewProb=0.5,simMethods=c("jaccard","geometric","invlogweighted","wt"),nThreads=1)

# Display summary information for simulation
print(re)

# Plot ROC curves for the result
plot(re)

# Get the AUC values for the result
re$AUC
}
\keyword{ROC; AUC}