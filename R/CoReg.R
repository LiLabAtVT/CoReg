# library("dynamicTreeCut")
# library("igraph")
# library("parallel")

#########################################################
# Main function for doing the network clustering
#########################################################

CoReg<-function(g,sim = "jaccard",minDegree = 1,nThreads=1){
  
  # Check input
  if(!is(g,"igraph")) stop("Argument g should be a graph object!")
  if(!is(nThreads,"numeric") || nThreads %% 1 != 0 || nThreads<1){ 
    stop("Argument nThreads should be a integer >=1")
  }
  if(!is(minDegree,"numeric") || minDegree %%1 != 0 || minDegree < 0){
    stop("minDegree should be an integer value >= 0")
  }
  
  # Get gene names
  g.degree<-degree(g)
  names(g.degree)<-as_ids(V(g))
  gene.names<-names(g.degree[g.degree>=minDegree])
  gene.idx<-which(g.degree>=minDegree)

  #gene.names<-V(g)$name
  
  # Choose between similarity measureament
  if(sim == "invlogweighted"){
    g.outsim.iw<-similarity.invlogweighted(g,gene.names,mode='out')
    g.outsim.iw<-g.outsim.iw[,gene.idx]
    g.insim.iw<-similarity.invlogweighted(g,gene.names,mode='in')
    g.insim.iw<-g.insim.iw[,gene.idx]
  }else if(sim == "jaccard"){
    g.outsim.iw<-similarity.jaccard(g,gene.names,mode='out')
    g.insim.iw<-similarity.jaccard(g,gene.names,mode='in')
  }else if(sim == "geometric"){
    g.outsim.iw<-getSimMatGeo(g,gene.names,mode='out',nThreads=nThreads)
    g.insim.iw<-getSimMatGeo(g,gene.names,mode='in',nThreads=nThreads)
  }else{
    stop("Please specify a correct similarity measure!")
  }
  
  # Add incoming and outgoing similarity together
  # Name the row and column for the similarity matrix
  g.sim.iw<-g.outsim.iw+g.insim.iw
  colnames(g.sim.iw)<-gene.names
  rownames(g.sim.iw)<-gene.names
  
  # Calculate a dissimilarity matrix based on similarity matrix 
  g.dsim.iw<- (max(g.sim.iw)-g.sim.iw)/max(g.sim.iw)
  diag(g.dsim.iw)<-0
  
  g.dist<-as.dist(g.dsim.iw)
  # Do hierarchical clustering
  re.hclust<-hclust(g.dist)
  
  # Cut the tree using dynamicTreeCut
  re.cut<-cutreeDynamic(re.hclust,distM = as.matrix(g.dist),method = "hybrid",minClusterSize=2,verbose = 0)
  
  # Sort the gene pairs based on the similarity score
  gene.names.mat<-t(combn(gene.names,2))
  sim<-g.sim.iw[gene.names.mat]
  rank.mat<-data.frame(gene.names.mat,sim,stringsAsFactors = F)
  rank.mat<-rank.mat[order(rank.mat[,3],decreasing=T),]
  colnames(rank.mat)<-c("gene1","gene2","similarity_score")
  
  # Generate output
  re.output<-list("module"=data.frame(ID=rownames(g.sim.iw),module=re.cut,stringsAsFactors = F),"similarity_matrix"=g.sim.iw,"rank"=rank.mat)
  
  class(re.output)<-"CoReg.result"
  return(re.output)
}

###############################################
# Some generic functions
###############################################
print.CoReg.result<-function(CoReg.result){
  numGene<-nrow(CoReg.result$module)
  numModule<-length(setdiff(unique(CoReg.result$module[,2]),0))
  cat(paste(numGene," genes in total\n",sep=""))
  cat(paste(numModule," module(s) were found"))
}

summary.CoReg.result<-function(CoReg.result){
  print(CoReg.result)
}

##############################################
# Calculate geometric index similarity matrix
# g: graph, mode: incoming or outgoing,
# n: number of threads, 
# vids: vertex ids for which similarity is
# calculated
##############################################
getSimMatGeo<-function(g,vids,mode,nThreads=1){
  gene.names<-vids
  if(mode == 'out'){
    nei.lst<-neighborhood(g,1,nodes=vids,mode='out')
  }else{
    nei.lst<-neighborhood(g,1,nodes=vids,mode="in")
  }
  
  # convert to normal vectors
  nei.lst<-sapply(nei.lst,as_ids)
  
  # remove node id itself
  nei.lst<-lapply(nei.lst,function(x){x<-x[-1]})
  
  # name the list
  names(nei.lst)<-gene.names
  nei.lst<<-nei.lst
  
  if(nThreads>1){
    # spec for parallelization
    cl<<-makeCluster(nThreads)
    clusterExport(cl,"nei.lst")
    
    # parallelized implementation for geometric index
    sim.mat<-parSapply(cl,
                       gene.names,
                       function(x){
                         sapply(
                           gene.names,
                           function(y){
                             length(intersect(nei.lst[[x]],nei.lst[[y]]))^2/(length(nei.lst[[x]])*length(nei.lst[[y]]))
                           }
                         )
                       }
    )
    stopCluster(cl)
  }else{
    
    # A non parallelized version if nThreads = 1
    sim.mat<-sapply(
      gene.names,
      function(x){
        sapply(
          gene.names,
          function(y){
            length(intersect(nei.lst[[x]],nei.lst[[y]]))^2/(length(nei.lst[[x]])*length(nei.lst[[y]]))
          }
        )
      }
    )
  }
  
  rm(nei.lst,envir=.GlobalEnv)
  # Set the cell without similarity score as zero
  sim.mat[is.na(sim.mat)]<-0
  return(sim.mat)
}