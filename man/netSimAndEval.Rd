\name{netSimAndEval}
\alias{netSimAndEval}
\alias{print.CoReg.NetSim}
\alias{summary.CoReg.NetSim}
\alias{plot.CoReg.NetSim}
\title{Generate simulated network and test}
\description{ This function performs generates simulated network(s) and use different module-finding methods to identify modules in the simulated networks. The performance of different methods will be evaluated by calculating NMI score between ground-truth modules and the algorithm identified modules}
\usage{
netSimAndEval(mSize,mNum,targetNum,auxNum,prob,nRep=1,testMethods=c("coregJac","coregInv","coregGeo","lp","wt","eb"),nThreads=1)
  \method{print}{CoReg.NetSim}()
  \method{summary}{CoReg.NetSim}()
  \method{plot}{CoReg.NetSim}()
}

\arguments{
  \item{mSize}{Size of module in simulated network(s)}
  \item{mNum}{Number of regulator gene in each module}
  \item{targetNum}{Number of target genes for each regualtor in the module}
  \item{auxNum}{Number of auxiliary genes in simulated network(s)}
  \item{prob}{Co-regulation probability. A vector of numeric values between 0 and 1. Larger value generates network with stronger co-regulation pattern}
  \item{nRep}{Number of replicates for each data point}
  \item{testMethods}{Methods used for finding modules. Their performance will be evaluated by NMI score. Available opotions are: "coregJac" (CoReg + Jaccard index), "coregGeo" (CoReg + geometric index),"coregInv" (CoReg + inverse log-weighted index),"lp" (label propagation),"wt" (walk trap),"eb" (edge betweenness)}
  \item{nThreads}{Number of threads for running the simulation. Only valid when "coregGeo" is selected in \code{Methods}}
}

\details{
  This function first calls \code{generateSimNet()} to generate simulated network(s) with pre-specified modular structure, and then runs different module-finding algorithms to identify modules. The correlation between pre-specified modules and algorithm identified modules is calculated using NMI score. The result can be plotted using generic function \code{plot()}}

\value{
  An object of class \code{CoReg.NetSim}
}

\seealso{
 \code{\link{generateSimNet}}
 \code{\link{rewSim}}
}
\author{Qi Song}
\references{Guimerà R, Sales-Pardo M, Amaral LAN: Module identification in bipartite and directed networks. Phys Rev E - Stat Nonlinear, Soft Matter Phys 2007, 76.
Danon L, Díaz-Guilera A, Duch J, Arenas A: Comparing community structure identification. J Stat Mech Theory Exp 2005, 2005:P09008–P09008.
}
\examples{
re<-netSimAndEval(10,5,20,100,0.5,testMethods=c("coregJac","lp","wt","eb"))

# See the summary
summary(re)

# Plot the evaluation result
plot(re)
}