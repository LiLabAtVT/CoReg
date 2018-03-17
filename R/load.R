#########################################################
# Load the data file and convert the directed network
# to a bipartite network.
#########################################################
# library(igraph)

networkFromFile<-function(file,sep = "",header = F,simple = T){
    el<-read.table(file=file,sep=sep,header=header,as.is=T,quote="")
    if(ncol(el)!=2){
      stop("Input file should contain two columns. Please check the format")
    }
    g<-graph.edgelist(as.matrix(el),directed = T)
    # Remove self loop and multiple edges
    if(simple == T){
        g<-simplify(g)
        return(g)
    }else if(!is(simple,"logical")){
        stop("Argument 'simple' should be 'True' or 'False'!")
    }
}

networkFromEdgeList<-function(edgeList,simple = T){
    if(ncol(edgeList)!=2 && !is.data.frame(edgeList)){
        stop("edgeList should be a data frame containing two columns. Please check the input")
    }
    g<-graph.edgelist(as.matrix(edgeList),directed = T)
    # Remove self loop and multiple edges
    if(simple == T){
        g<-simplify(g)
        return(g)
    }else if(!is(simple,"logical")){
        stop("Argument 'simple' should be 'True' or 'False'!")
    }
}
#########################################################
# Convert a directed network to a bipartite one
# Note this is actually stored as undirected igraph
# graph object.
#########################################################
directedToBipartite<-function(g){
    if(!is(g,"igraph")) stop("Argument g should be an igraph object")
    el<-as_edgelist(g)
    el[,1]<-paste(el[,1],'_h',sep='')
    el[,2]<-paste(el[,2],'_t',sep='')
    g<-graph.edgelist(as.matrix(el),directed = F)
    return(g)  
}