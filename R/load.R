#########################################################
# Load the data file and convert the directed network
# to a bipartite network.
#########################################################
library(igraph)

loadNetwork<-function(file,sep = "",header = F,mode = 'd',simple = T){
    el<-read.table(file=file,sep=sep,header=header,as.is=T,quote="")
    if(mode == 'b'){
        g<-graph.edgelist(as.matrix(el),directed = F)
        g<-directedToBipartite(g)
    }else if(mode == 'd'){
        g<-graph.edgelist(as.matrix(el),directed = T)  
    }else{
        stop("Please specify correct mode!")
    }
    
    # Remove self loop and multiple edges
    if(simple == T){
        g<-simplify(g)
        return(g)
    }else if(!is(simple,"logical")){
        stop("Argument 'simple' should be logical!")
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