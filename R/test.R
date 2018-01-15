# library("igraph")

####################################################################################
# First convert directed to bipartite and then run community detection algorithms on
# top of if
#####################################################################################
test<-function(g,test.method){
    bipartite
    switch(test.method,
             "lp"={
               g.bi<-directedToBipartite(g)
               cluster<-label.propagation.community(g.bi) 
             },
             "wt"={
               g.bi<-directedToBipartite(g)
               cluster<-walktrap.community(g.bi,steps=2)
             },
             "eb"={
               g.bi<-directedToBipartite(g)
               cluster<-edge.betweenness.community(g.bi)
             },
             "coregJac"={
               cluster<-CoReg(g,"jaccard",minDegree = 0)$module
             },
             "coregInv"={
               cluster<-CoReg(g,"invlogweighted",minDegree = 0)$module
             },
             "coregGeo"={
               cluster<-CoReg(g,"geometric",minDegree = 0)$module
             }
           )
    
    if(exists("g.bi",envir = )){
      re.mat<-getcom(cluster$membership,g.bi)
      re.mat<-flatModule(re.mat)
      return(re.mat)
    }else{
      return(cluster)
    }
}

###############################################################
# Function to get the module assignment for lp,wt and eb
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
        if(length(k)==1){outmat[each,2]<-k}
        else{outmat[each,2]<-paste(k,collapse=';')}
    }
    return(outmat)
}

#######################################################
# Some gene might have multiple module assignments
# This function split the multiple module assignments
# of the same gene into different row
#######################################################
flatModule<-function(module){
  result<-data.frame("ID"=NA,"module"=NA)
  for(i in 1:nrow(module)){
    module_list<-as.numeric(unlist(strsplit(module[i,2],";")))
    gene<-rep(module[i,1],times=length(module_list))
    new<-data.frame("ID" = gene,"module" = module_list)
    result<-rbind(result,new)
  }
  result<-result[-1,]
  return(result)
}
