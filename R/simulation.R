#########################################################
# This modules consists of functions that do network global
# rewiring simulation and caculation of similarities 
# between original graph and the duplicated one. 
#########################################################
# library(igraph)
# library(ggplot2)

########################################################
# This is the main function of this module
########################################################
rewSim<-function(g,nDup,dDup,rewProb,methods=c(),nRep = 1,nThreads =1){
    
    ########################################################
    # Check input
    ########################################################
    if(!is(g,"igraph")) stop("g argument only takes a graph object!")
    if(!is(nDup,"numeric") || nDup %% 1 != 0 || nDup<1 || nDup>gorder(g)){
      stop("Argument nDup should be a numeric value between 1 and number of nodes in the graph!")
    }
    if(!is(nRep,"numeric") || nRep %% 1 != 0 ||nRep<1){
      stop("Argument nRep should be a positive integer >= 1!")
    }
    if(!is(nThreads,"numeric") || nThreads %% 1 != 0 || nThreads<1){
      stop("Argument nThreads should be an integer >= 1!")
    }
    if(!is(rewProb,"numeric") || TRUE %in% rewProb>1 || TRUE %in% rewProb<0){ 
      stop("Argument rewProb should be a numeric vector containing values between 0 and 1!")
    }
    allMethods = c("lp","wt","eb","coregInv","coregJac","coregGeo")  
    for(method in methods){
      if (!(method %in% allMethods)){
        stop(paste("'",method,"'"," is not a valid method name",sep=""))  
      }
    }
    
    ########################################################
    # Part I.Duplicate and rewire the graph then partition
    # the network using our clustering algorithm
    ########################################################

    # Duplicate the graph  
    re<-networkDup(g,nDup,dDup)
    g.dup<-re$graph
    actualDup<-re$nDup
    
    # Create scores matrix/matrices
    scoreLst<-list()
    
    for(method in methods){
        scoreLst[[paste(method,"score",sep = "-")]]=matrix(nrow=length(rewProb),row.names<-rewProb,ncol=nRep)
    }
  
    # How many reps you want for each data point?
    for(iter in 1:nRep){
        # Rewire the duplicates, get the score for the merged graph
        for(i in 1:length(rewProb)){
            g.dup.rewired<-networkRew(g.dup,rate=rewProb[i])
            g.merge<-g+g.dup.rewired
    
            ########################################################
            # Part II. Do exactly the same thing as part I but use
            # lp, wt and eb clustering algorithms instead
            ########################################################
          
            for(method in methods){
                scoreLst[[paste(method,"score",sep="-")]][i,iter]<-RRS_Eval(test(g.merge,test.method=method))
            }
        }
    }
    
    # reformat score.lst replace the value with its mean and standard deviation
    scoreLst<-reformatResult(scoreLst,nRep = nRep)
     
    # Create results matrix
    results.output<-do.call(data.frame,scoreLst)
    
    # Format column names
    colnames(results.output)<-gsub("\\.","-score-",colnames(results.output))
    
    # Make the data frame for plotting purpose
    probColumn<-rep(rewProb, times=length(scoreLst))
    methodColumn<-rep(names(scoreLst), each = length(rewProb))
    results.plot<-data.frame(rewProb = probColumn, do.call(rbind,scoreLst), method = methodColumn)
    
    re<-list(plotData=results.plot,evalResult=results.output,nDup=actualDup,dDup=dDup,methods=methods,nRep=nRep,rewProb=rewProb)
    
    class(re)<-"CoReg.RewSim"
    return(re)
}

##########################################################################
# Plot method for CoRegRew object. CoRegRew object
# is returned by rewSim. This function can plot performance of selected 
# method(s), evaluated by RRS score, which is computed during the rewiring
# simulation on real networks.
############################################################################
plot.CoReg.RewSim<-function(CoReg.RewSim){
  figure<-ggplot(CoReg.RewSim$plotData,aes(x=rewProb,y=mean,colour=method)) +
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.01) +
    geom_line() +
    geom_point() +
    xlab("Rewiring Probability") +
    ylab("Rewiring recall score") +
    ylim(0,1) +
    theme(text=element_text(size=16),axis.text=element_text(size=16),axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))
  print(figure)
}

print.CoReg.RewSim<-function(CoReg.RewSim){
  cat("CoReg.RewSim object\n")
  cat("-------------Summary of rewiring simulation-------------\n")
  cat(CoReg.RewSim$nDup,"nodes with at least",CoReg.RewSim$dDup,"neighbors were duplicated.\n")
  cat("Each data point is calculated as average of",CoReg.RewSim$nRep,"run(s).\n")
  usedMethods<-do.call(paste,c(as.list(CoReg.RewSim$methods),sep=" "))
  cat("Clustering methods used:",usedMethods)
  cat(".\n")
  rewProb = do.call(paste,c(as.list(CoReg.RewSim$rewProb),sep=" "))
  cat("Rewiring probability tested:",rewProb)
  cat(".\n")
}

summary.CoReg.RewSim<-function(CoReg.RewSim){
  print(CoReg.RewSim)
}

#########################################################
# Function for duplicating original graph, the duplicates
# will be returned
#########################################################
networkDup<-function (g,nDup,dDup) {
  vertex.lst <- names(which(degree(g, V(g)$name, mode = "all") >= dDup))
  
  if (nDup > length(vertex.lst)) {
    duplicates.names <- vertex.lst
  }
  else {
    duplicates.names <- sample(vertex.lst, size = nDup, replace = F)
  }
  g.sim <- make_empty_graph(n = 0, directed = T)
  for (node.name in duplicates.names) {
    nei.outnames <- neighbors(g, v = node.name, mode = "out")$name
    nei.innames <- neighbors(g, v = node.name, mode = "in")$name
    node.name <- paste("dup-", node.name, sep = "")
    in.edges.lst <- c(t(cbind(nei.innames, rep(node.name, 
                                               times = length(nei.innames)))))
    out.edges.lst <- c(t(cbind(rep(node.name, times = length(nei.outnames)), 
                               nei.outnames)))
    nei.innames <- setdiff(nei.innames, V(g.sim)$name)
    nei.outnames <- setdiff(nei.outnames, V(g.sim)$name)
    g.sim <- add.vertices(g.sim, nv = (1 + length(c(nei.innames, nei.outnames))), name = c(node.name, nei.innames,nei.outnames))
    g.sim <- add.edges(g.sim, edges = c(in.edges.lst, out.edges.lst))
  }
  return(list(graph=g.sim,nDup=length(duplicates.names)))
}

#########################################################
# Evaluate the clustering results. Take out the duplicates
# and see whether they stay in the same module with the
# original nodes
#########################################################
RRS_Eval<-function(modules){
    
    # n1 is the total number of genes and n2 is the number of duplicates genes 
    n1<-length(modules[,1])
    n2<-length(grep("^dup-",modules[,1]))
  
    score.lst<-c()
    
    # Penalty coefficient for each cluster
    penalty.lst<-c()
    
    # Count number of duplicated nodes that stay in the same cluster with the original nodes
    # and calculate the penalty coefficient for each node
    for(i in unique(modules[,2])){
        
        geneset<-as.character(modules[which(modules[,2]==i),1])
        duplicates<-geneset[grep("^dup-",geneset)]
        
        if(length(duplicates) > 0){
            duplicates<-unlist(strsplit(duplicates,"dup-"))
            duplicates<-duplicates[-which(duplicates=="")]
            score.lst<-c(score.lst,rep(1,times=length(intersect(duplicates,geneset))),rep(0,times=length(setdiff(duplicates,geneset))))
            penalty.lst<-c(penalty.lst,rep(n1/length(geneset),times=length(duplicates)))
        }
    }
  
    # Add penalty coefficient to the score
    score.lst<-score.lst*penalty.lst
    score<-sum(score.lst)/(n1*n2/2)
    
    return(score)
}

#########################################################
# Function for rewiring the network
#########################################################
networkRew<-function(g,rate){
    
    #####################################################
    # Part I. Rewire the directed network
    #####################################################
    if(is.directed(g)){
        all.names<-V(g)$name
        duplicates.names<-all.names[grep("^dup-",all.names)]
        el.duplicates.out<-matrix(ncol=2);el.duplicates.in<-matrix(ncol=2)
        
        # Generate an outgoing edgelis and an incoming edgelist
        # Get all the outgoing neighbors and incoming neighbors
        # Then create the outgoing and incoming edge lists
        for(gene.name in duplicates.names){
            duplicates.outnames<-neighbors(g,gene.name,mode = "out")$name  
            duplicates.innames<-neighbors(g,gene.name,mode = "in")$name
            if(length(duplicates.outnames)>0){
                el.new.out<-matrix(c(rep(gene.name,times=length(duplicates.outnames)),duplicates.outnames),ncol=2)
                el.duplicates.out<-rbind(el.duplicates.out,el.new.out)
            }
            if(length(duplicates.innames)>0){
                el.new.in<-matrix(c(duplicates.innames,rep(gene.name,times=length(duplicates.innames))),ncol=2)
                el.duplicates.in<-rbind(el.duplicates.in,el.new.in)
            }
        }
        
        # Remove first empty row
        el.duplicates.out<-el.duplicates.out[-1,]
        el.duplicates.in<-el.duplicates.in[-1,]
        
        # Exchange edges for incoming edge list and outgoing edge list
        lst.out<-el.duplicates.out[,1]
        lst.in<-el.duplicates.in[,2]
    
        for(i in 1:length(lst.out)){
            if(runif(1,0,1) <= rate && length(which(lst.out != lst.out[i])) > 0){
                index<-sample(which(lst.out != lst.out[i]),size=1)
                tmp<-el.duplicates.out[i,2]
                el.duplicates.out[i,2]<-el.duplicates.out[index,2]
                el.duplicates.out[index,2]<-tmp
            }
        }
    
        
        for(i in length(lst.in)){
            if(runif(1,0,1) <= rate && length(which(lst.in != lst.in[i])) > 0){
                index<-sample(which(lst.in != lst.in[i]),size=1)
                tmp<-el.duplicates.in[i,2]
                el.duplicates.in[i,2]<-el.duplicates.in[index,2]
                el.duplicates.in[index,2]<-tmp
            }
        }
    
        # Create a graph using the two edge lists and return it
        g.rewired<-make_empty_graph(n=0,directed=T)
        g.rewired<-add.vertices(g.rewired,nv=length(all.names),name = all.names)
        seq.duplicates.in<-c(t(el.duplicates.in))
        seq.duplicates.out<-c(t(el.duplicates.out))
        g.rewired<-add.edges(g.rewired,edges = c(seq.duplicates.in,seq.duplicates.out))
        return(g.rewired)    
    }else{
        
        #######################################################################
        # Part II. Rewire the undirected (bipartite) network
        #######################################################################
        all.names<-V(g)$name
        duplicates.names<-all.names[grep("^dup-",all.names)]
        el.duplicates<-matrix(ncol=2)
        
        for(gene.name in duplicates.names){
            duplicates.neighbors<-neighbors(g,gene.name,mode = "all")$name  
            el.new<-matrix(c(rep(gene.name,times=length(duplicates.neighbors)),duplicates.neighbors),ncol=2)
            el.duplicates<-rbind(el.duplicates,el.new)
        }
        
        # Remove first empty row
        el.duplicates<-el.duplicates[-1,]
        el.duplicates.original<-el.duplicates
        
        # Exchange edges
        if(!is.matrix(el.duplicates)) el.duplicates<-matrix(el.duplicates,ncol=2)
        lst.genes<-el.duplicates[,1]
        for(i in 1:length(lst.genes)){
            if(runif(1,0,1) <= rate && length(which(lst.genes != lst.genes[i])) > 0){
                index<-sample(which(lst.genes != lst.genes[i]),size=1)
                tmp<-el.duplicates[i,2]
                el.duplicates[i,2]<-el.duplicates[index,2]
                el.duplicates[index,2]<-tmp
            }
        }
        
        # Create a graph using the edge list and return it
        g.rewired<-make_empty_graph(n=0,directed = F)
        g.rewired<-add.vertices(g.rewired,nv=length(all.names),name = all.names)
        seq.duplicates<-c(t(el.duplicates))
        g.rewired<-add.edges(g.rewired,edges = c(seq.duplicates))

        return(g.rewired)    
    }
}

##############################################################
# A wrapper function that does the following things in order
#     1. Generate a simulated network (including the true module
#        partition)
#     2. Run different module algorithm on the simulated network
#     3. Compare the algorithm-identified modules to the true module
#        partition and calculates the score
##############################################################
netSimAndEval<-function(mSize,mNum,targetNum,auxNum,prob,nRep=1,testMethods=c("coregJac","coregInv","coregGeo","lp","wt","eb"),nThreads=1){
  
  # Check the input
  if(!is(mSize,"numeric") || mSize %%1 != 0 || mSize <= 0){
    stop("mSize should be a integer >= 1")
  }
  if(!is(mNum,"numeric") || mNum %%1 != 0 || mNum <= 0){
    stop("mNum should be a integer >= 1")
  }
  if(!is(targetNum,"numeric") || targetNum %%1 != 0 || targetNum <= 0){
    stop("targetNum should be a integer >= 1")
  }
  if(!is(auxNum,"numeric") || auxNum %%1 != 0 || auxNum <= 0){
    stop("auxNum should be a integer >= 1")
  }
  if(!is(prob,"numeric") || TRUE %in% prob>1 || TRUE %in% prob<0){
    stop("prob should only contain numeric between 0 and 1")
  }
  if(!is(nRep,"numeric") || nRep %%1 != 0 || nRep <= 0){
    stop("nRep should be an integer >= 1")
  }
  
  # Create scores matrix/matrices
  scoreLst <- list()
  for(method in unique(testMethods)){
    if(method !="lp" && method != "wt" && method != "eb" && method != "coregJac" && method != "coregInv" && method != "coregGeo" ){
      stop("Please specify correct method name for compareMethods!")  
    }
    scoreLst[[method]] = matrix(nrow=length(prob),row.names<-prob,ncol=nRep)
  }
  
  # Test different algorithm using different prob and there are nRep replicates
  # for each prob
  for(i in 1:length(prob)){
    for(j in 1:nRep){
      
      # Genereate a simulated network with ground truth known (truePartition)
      simNet<-generateSimNet(mSize,mNum,targetNum,auxNum,prob[i])
      g<-graph.edgelist(apply(simNet$el,2,as.character))
      
      # Evaluate the module finding result only for the first mSize*mNum genes
      geneToEvaluate<-as.character(1:(mSize*mNum))
      truePartition<-simNet$modulePartition
      
      # Run different algorithm. Compare the result to the truePartition
      for(method in unique(testMethods)){
        geneAndPartition<-test(g,method)
        
        # To make the length of module assignment equal to each other
        # Only one module assignment is selected for each gene
        scoreLst[[method]][i,j]<-compare(geneAndPartition[match(geneToEvaluate,geneAndPartition[,"ID"]),"module"],truePartition,method="nmi")
      }
    }
  }
  scoreLst<-reformatResult(scoreLst,nRep)
  
  # Create results matrix
  results.output<-do.call(data.frame,scoreLst)
  
  # Format column names
  colnames(results.output)<-gsub("\\.","-score-",colnames(results.output))
  
  # Make the data frame for plotting purpose
  probColumn<-rep(prob, times=length(scoreLst))
  methodColumn<-rep(names(scoreLst), each = length(prob))
  results.plot<-data.frame(prob = probColumn, do.call(rbind,scoreLst), method = methodColumn)
  
  re<-list(
            plotData=results.plot,evalResult=results.output,mSize=mSize,
            mNum=mNum,targetNum=targetNum,auxNum=auxNum,prob=prob,methods=testMethods,
            nRep=5
          )
  class(re)<-"CoReg.NetSim"
  return(re)
}

##########################################################################
# Plot method for CoReg.NetSim object. CoReg.NetSim object
# is returned by netSimAndEval. This function can
# plot the clustering performance of selected method(s)
# , evaluated by NMI score between module assignment from
# the algorithm and pre-specified module assignment during
# simulation
############################################################################
plot.CoReg.NetSim<-function(CoReg.NetSim){
  
  figure<-ggplot(CoReg.NetSim$plotData,aes(x=prob,y=mean,colour=method)) +
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.01) +
    geom_line() +
    geom_point() +
    xlab("Probability") +
    ylab("NMI") +
    ylim(0,1) +
    theme(text=element_text(size=16),axis.text=element_text(size=16),axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))
  print(figure)
}

print.CoReg.NetSim<-function(CoReg.NetSim){
  cat("CoReg.NetSim object\n")
  cat("-------------Summary of simulated networks-------------\n")
  # rewProb = unique(CoReg.RewSim$evalResult)
  cat(CoReg.NetSim$mNum,"modules of size",CoReg.NetSim$mSize,"were generated.\n")
  cat("Each regulator has",CoReg.NetSim$targetNum,"targets.\n")
  cat("Each data point is calculated as average of",CoReg.NetSim$nRep,"run(s).\n")
  usedMethods<-do.call(paste,c(as.list(CoReg.NetSim$methods),sep=" "))
  cat("Clustering methods used:",usedMethods)
  cat(".\n")
  prob = do.call(paste,c(as.list(CoReg.NetSim$prob),sep=" "))
  cat("co-regulation probability tested:",prob)
  cat(".\n")
}
summary.CoReg.NetSim<-function(CoReg.NetSim){
  print(CoReg.NetSim)
}

##############################################################
# Function for generating a simulated network with ground
# truth known (module assignment). The module assignment of
# the first mSize*mNum genes are returned
##############################################################
generateSimNet<-function(mSize,mNum,targetNum,auxNum,prob){
  
  # Check the input
  if(!is(mSize,"numeric") || mSize %%1 != 0 || mSize <= 0){
    stop("mSize should be a integer >= 1")
  }
  if(!is(mNum,"numeric") || mNum %%1 != 0 || mNum <= 0){
    stop("mNum should be a integer >= 1")
  }
  if(!is(targetNum,"numeric") || targetNum %%1 != 0 || targetNum <= 0){
    stop("targetNum should be a integer >= 1")
  }
  if(!is(auxNum,"numeric") || auxNum %%1 != 0 || auxNum <= 0){
    stop("auxNum should be a integer >= 1")
  }
  if(!is(prob,"numeric") || TRUE %in% prob>1 || TRUE %in% prob<0){
    stop("prob should only contain numeric value from 0 to 1")
  }
  
  # Generate the node ID and module ID for each node
  allNodeID<-seq(mSize*mNum+auxNum)
  
  # Only the module assignments of the first mSize*mNum genes
  # are return
  modulePartition<-rep(seq(mNum),each = mSize)
  
  # Select the pool of candidates for each module
  allCandidatePool<-list()
  for(moduleID in unique(modulePartition)){
      # currentModuleNodeID
      allCandidatePool[[moduleID]]<-sample(allNodeID,targetNum)
  } 
  
  # A vector for storing the edgelist of the simulated network
  el<-c()
  
  # For each regulator, either select a target from the candidate pool (with 
  # probability specified by prob), or select that from out side the pool.
  # mNum+1 is the module ID for auxilary nodes
  for(nodeID in 1:length(modulePartition)){
    
      # Decide how many genes should be selected from candidate pool
      # And how many should be selected from outside the candidate pool
      decisionMaker<-runif(targetNum)
      NumberInPool = length(decisionMaker[decisionMaker<=prob])
      NumberOutPool = targetNum - NumberInPool
      
      # Rmove current node from the pool to avoid self-interaction
      candidateInPool = setdiff(allCandidatePool[[modulePartition[nodeID]]],nodeID)
      candidateOutPool = setdiff(allNodeID,c(candidateInPool,nodeID))
      
      # After removing the current node itself. The number of candidates in the pool
      # might be less than the desirable number of targets. In that case, we 
      # select all the genes in the pool
      NumberInPool = ifelse(length(candidateInPool)>NumberInPool,NumberInPool,length(candidateInPool))
      NumberOutPool = ifelse(length(candidateOutPool)>NumberOutPool,NumberOutPool,length(candidateOutPool))
      
      # select the targets from the pool
      targetInPool<-sample(candidateInPool,NumberInPool)
      targetOutPool<-sample(candidateOutPool,NumberOutPool)
      el<-c(el,as.vector(matrix(c(rep(nodeID,NumberInPool+NumberOutPool),targetInPool,targetOutPool),nrow=2,byrow=T)))
  }
  el<-matrix(el,ncol=2,byrow=T)

  return(list("el"=el,"modulePartition"=modulePartition))
}

######################################################
# A helper function for netSimAndEval and rewSim
# This function reformats the result from
# the simulation to have clean output
######################################################
reformatResult<-function(scoreLst,nRep){
  
  # reformat scoreLst replace the value with its mean and standard deviation
  for(method in names(scoreLst)){
    score.mean<-apply(scoreLst[[method]],1,mean)
    if(nRep<=1){
      score.sd<-rep(0,times=nrow(scoreLst[[method]])) 
    }else{
      score.sd<-apply(scoreLst[[method]],1,sd)
    }
    scoreLst[[method]]=matrix(ncol=2,c(score.mean,score.sd))
    colnames(scoreLst[[method]])<-c("mean","sd")
  }
  
  return(scoreLst)
}