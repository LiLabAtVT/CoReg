\name{CoReg}
\alias{CoReg}
\alias{print.CoReg.result}
\alias{summary.CoReg.result}
\title{Identify Co-regulators in transcription gene network}
\description{This is the main function for identifying the co-regulators in transcription gene network}
\usage{
CoReg(g, sim = "jaccard", minDegree = 1, nThreads = 1)
\method{print}{CoReg.result}(CoReg.result)
\method{summary}{CoReg.result}(CoReg.result)
}
\arguments{
\item{g}{ The input transcription network. This should be an "igraph" object returned by the function \code{"networkFromFile"} or \code{"networkFromEdgeList"}  }
\item{sim}{ Similarity index for computing the similarity between genes. Available options are: "jaccard", "geometric", "invlogweighted". Default value is "jaccard". }
\item{minDegree}{ The minimum degree for genes that need similarity calculation. Genes with degree smaller than minDegree will not be computated and therefore will not be considered as co-regulator to any other genes}
\item{nThreads}{ Number of threads to be used for similarity computation. Only valid when \code{sim = "geometric"}. }
\item{CoReg.result}{ \code{CoReg.result} object}
}
\details{
  This function takes an "igraph" object produced by \code{"networkFromFile()"} or \code{"networkFromEdgeList()"}. The returned object includes the similarity matrix, rank and co-regulatory module assignment computed using the similarity index specified in \code{sim} argument. CoReg first calculates the pairwise similarity among all the genes, then performs hierarchical clustering and dynamic tree cut to obtain the final module assignment. See the following paper for more details about dynamic tree cut algorithm: Langfelder P, Zhang B, Horvath S (2007) Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. Bioinformatics 2008 24(5):719-720.
}
\value{
An object of class \code{CoReg.result}, a list including elements
\item{module}{ a dataframe reporting the module ID for each gene }
\item{similarity_matrix}{ the similarity matrix for all the input genes  }
\item{rank}{ the rank of similarity for all the gene pairs}
}
\seealso{
  \code{\link{networkFromFile}}
  \code{\link{networkFromEdgeList}}
}
\author{Qi Song}
\references{Langfelder P, Zhang B, Horvath S (2007) Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. Bioinformatics 2008 24(5):719-720.}
\examples{
data(athNet)
re <- CoReg(athNet)

# Display summary information for the identfied co-regulatory modules
summary(re)
print(re)

# Get the similarity matrix from the result
re$similarity_matrix

# Get the simialrity rank from the result
re$rank

# Get the module assignment from the result
re$module
}
\keyword{co-regulator}