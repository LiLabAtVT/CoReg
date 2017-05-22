\name{loadNetwork}
\alias{loadNetwork}
\title{Load transcription network}
\description{Function for loading a transcription network. Input should be an edge list}
\usage{
loadNetwork(file,sep = "",header = F,mode = 'd',simple = T)
}
\arguments{
\item{file}{ the file name for the network edge list}
\item{sep}{ the delimiter separating the columns in the input file }
\item{header}{ a boolean value indicating whether the input file contains header. \code{TRUE} for yes, \code{FALSE} for no. Default value is \code{FALSE} }
\item{mode}{ whether input network is a directed or undirected network. 'd' for directed network, 'u' for undirected network.}
\item{simple}{ a boolean value indicating whether to load the network as simple graph. \code{TRUE} for yes, \code{FALSE} for no. Default value is \code{TRUE}. \code{simple = TRUE} means the self-loops and multiple edges will be removed once the network is loaded. }
}
\details{
  This function loads the network from file. Network should be fomated into a two-column edge list in the file. The first column represents the gene ID for transcription factor and the second column represents the gene ID for the correspoinding target. Two columns should be separated by a delimiter (tab, comma, space, etc.)
}
\value{
An "igraph" object storing the network
}
\seealso{
  \code{\link{CoReg}}
}
\author{Qi Song}
\references{ Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. }
\examples{
}
\keyword{load network}