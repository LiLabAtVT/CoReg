library("igraph")

#########################################################################
# Convert directed to bipartite and run community detection algorithms on
# top of if
#########################################################################
test<-function(g,test.method){
    
    g.bi<-directedToBipartite(g)
    
    if(test.method=="lp")    cluster<-label.propagation.community(g.bi)
    if(test.method=="wt")    cluster<-walktrap.community(g.bi)
    if(test.method=="eb")    cluster<-edge.betweenness.community(g.bi)
    
    re.mat<-getcom(cluster$membership,g.bi)
    rownames(re.mat)<-c()
    return(re.mat)
}

###############################################################
# Function to get the original names back
###############################################################
getcom<-function(test,g){
    # test is the membership assignment
    names(test)<-V(g)$name
    test.mat<-cbind(test,names(test)) # two colume matrix, 
    # col 1: membership assignment
    # col 2: edge name
  
    # get the original gene names back
    getname<-function(x){strsplit(x,'_')[[1]][1]}
    allname<-unique(unlist(lapply(test.mat[,2],getname)))
  
    # include original name into the matrix
    test.mat1<-cbind(test.mat,lapply(test.mat[,2],getname))
  
    outmat<-data.frame(name=allname,group="",stringsAsFactors = F)
    rownames(outmat)<-outmat[,1]
  
    # prepare output matrix, with gene name and gene module assignment
    # for gene belong to two modules, assign both module names to the same gene
    for(each in allname)
    {
        k<-unname(unlist(test.mat1[test.mat1[,3]==each,1]))
        k<-sort(unique(k))
        #print(k);print(length(k));print(paste(k,collapse=';'))
        if(length(k)==1){outmat[each,2]<-k}
        else{outmat[each,2]<-paste(k,collapse=';')}
    }
    return(outmat)
}