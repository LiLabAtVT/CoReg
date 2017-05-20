#########################################################
# This modules consists of functions that do network global
# rewiring simulation and caculation of similarities 
# between original graph and the duplicated one. 
#########################################################
library(igraph)
library(ggplot2)
#source("test.R")
#source("clustering.R")
#source("load.R")

########################################################
# This is the main function of this module
########################################################
networkSim<-function(g,nDup,dDup,rewProb,compareMethods=c(),nRep = 1,nThreads =1){
    
    ########################################################
    # Check input
    ########################################################
    if(!is(g,"igraph")) stop("g argument only takes a graph object!")
    if(!is(nDup,"numeric") || nDup<1 || nDup>gorder(g)) stop("Argument nDup should be a numeric value from 1 to number of nodes in the graph!")
    if(!is(nRep,"numeric") || nRep<1) stop("Argument nRep should be an integer value >= 1!")
    if(!is(nThreads,"numeric") || nThreads<1) stop("Argument nThreads should be an integer value >= 1!")
    
    ########################################################
    # Part I.Duplicate and rewire the graph then partition
    # the network using our clustering algorithm
    ########################################################

    # Duplicate the graph  
    g.dup<-networkDup(g,nDup,dDup)
    
    # Create scores matrix/matrices
    score.lst<-list("coreg.inv.score"=matrix(nrow=length(rewProb),row.names<-rewProb,ncol=nRep),"coreg.jac.score"=matrix(nrow=length(rewProb),row.names<-rewProb,ncol=nRep),"coreg.geo.score"=matrix(nrow=length(rewProb),row.names<-rewProb,ncol=nRep))
    
    for(method in compareMethods){
        if(method !="lp" && method != "wt" && method != "eb" ){
            stop("Please specify correct method name for compareMethods!")  
        }
        score.lst[[paste(method,"score",sep = ".")]]=matrix(nrow=length(rewProb),row.names<-rewProb,ncol=nRep)
    }
  
    # How many reps you want for each data point?
    for(iter in 1:nRep){
        
        # Rewire the duplicates, get the score for the merged graph
        for(i in 1:length(rewProb)){
            if(!is(rewProb[i],"numeric") || rewProb[i]>1 || rewProb[i]<0) stop("Argument nDup should be a numeric value from 0 to 1!")
            g.dup.rewired<-networkRew(g.dup,rate=rewProb[i])
            g.merge<-g+g.dup.rewired
    
            ########################################################
            # Part II. Do exactly the same thing as part I but use
            # lp, wt and eb clustering algorithms instead
            ########################################################
            
            # Convert the merged graph to a bipartite graph so that the three algorithms could use
            g.merge.bi<-directedToBipartite(g.merge)
            module.test.rewired<-list()
            module.test.rewired[["coreg.inv"]]<-CoReg(g.merge,"invlogweighted")$module
            module.test.rewired[["coreg.jac"]]<-CoReg(g.merge,"jaccard")$module
            module.test.rewired[["coreg.geo"]]<-CoReg(g.merge,"geometric",nThreads=nThreads)$module
            score.lst[["coreg.inv.score"]][i,iter]<-networkEval(module.test.rewired[["coreg.inv"]])
            score.lst[["coreg.jac.score"]][i,iter]<-networkEval(module.test.rewired[["coreg.jac"]])
            score.lst[["coreg.geo.score"]][i,iter]<-networkEval(module.test.rewired[["coreg.geo"]])
            
            for(method in compareMethods){
                module.test.rewired[[method]]<-test(g.merge.bi,test.method=method)  
                module.test.rewired[[method]]<-get_module(module.test.rewired[[method]])
                score.lst[[paste(method,"score",sep=".")]][i,iter]<-networkEval(module.test.rewired[[method]])
            }
        }
    }
    
    # reformat score.lst replace the value with its mean and standard deviation
    for(method in names(score.lst)){
        score.mean<-apply(score.lst[[method]],1,mean)
        if(nRep<=1){
            score.sd<-rep(0,times=nrow(score.lst[[method]])) 
        }else{
            score.sd<-apply(score.lst[[method]],1,sd)
        }
        score.lst[[method]]=matrix(ncol=2,c(score.mean,score.sd))
        colnames(score.lst[[method]])<-c(paste(method,"mean",sep="."),paste(method,"sd",sep="."))
    }
        
    # Create results matrix
    results.output<-data.frame(rewProb=rewProb,score.lst[["coreg.inv.score"]],score.lst[["coreg.jac.score"]],score.lst[["coreg.geo.score"]])
    results.plot<-data.frame(rew=rewProb,mean=score.lst[["coreg.inv.score"]][,1],sd=score.lst[["coreg.inv.score"]][,2],method=rep("CoReg_inv",times=length(rewProb)))
    results.plot<-rbind(results.plot,data.frame(rew=rewProb,mean=score.lst[["coreg.jac.score"]][,1],sd=score.lst[["coreg.jac.score"]][,2],method=rep("CoReg_jaccard",times=length(rewProb))))
    results.plot<-rbind(results.plot,data.frame(rew=rewProb,mean=score.lst[["coreg.geo.score"]][,1],sd=score.lst[["coreg.geo.score"]][,2],method=rep("CoReg_geometric",times=length(rewProb))))
    
    for(method in names(score.lst)){
        if(method!="coreg.inv.score" && method!="coreg.jac.score" && method!="coreg.geo.score" ){
            results.output<-data.frame(results.output,score.lst[[method]])   
            if(method == "lp.score") method.name="label propagation"
            if(method == "eb.score") method.name="edge betweenness"
            if(method == "wt.score") method.name="walk trap"
            results.plot<-rbind(results.plot,data.frame(rew=rewProb,mean=score.lst[[method]][,1],sd=score.lst[[method]][,2],method=rep(method.name,times=length(rewProb))))
        }
    }
    
    re<-list(plotData=results.plot,evalResult=results.output)
    class(re)<-"CoReg.rewSim"
    return(re)
}

#######################################################
# Plot simulation results. If the parameter 'test'
# is enabled, then four curves of scores of each 
# algorithm will display. Otherwise only one curve will 
# show up.
########################################################
plot.CoReg.rewSim<-function(rewSim){
    
    figure<-ggplot(rewSim$plotData,aes(x=rew,y=mean,colour=method)) +
               geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.01) +
               geom_line() +
               geom_point() +
               xlab("Rewiring probability") +
               ylab("Scores") +
               ylim(0,1) +
               theme(text=element_text(size=16),axis.text=element_text(size=16),axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))
    print(figure)
}

print.CoReg.rewSim<-function(rewSim){
    print("CoReg.rewSim object")
}

getRewScore<-function(rewSim){
    if(!is(rewSim,"CoReg.rewSim")) stop("Argument rewSim should be an CoReg.rewSim object!")
    rewSim$evalResult
}

#########################################################
# Function for duplicating original graph, the duplicates
# will be returned
#########################################################
networkDup<-function(g,nDup,dDup){
    
    # Get genes whose degree are greater than dDup
    vertex.lst<-names(which(degree(g,V(g)$name,mode = "all")>=dDup))
    
    # Sample nDup nodes. If number of nodes are less than nDup, all the nodes will be taken
    if(nDup>length(vertex.lst)){
        duplicates.names<-vertex.lst    
    }else{
        duplicates.names<-sample(vertex.lst,size = nDup,replace = F)
    }
    
    # Get all the node names for the duplicated subgraph
    allNodeName<-unique(unlist(lapply(neighborhood(g,1,duplicates.names,mode="all"),as_ids)))
    
    # Generate duplicated subgraph
    g.sim<-induced_subgraph(g,allNodeName,"auto")
    
    # Replace the node names to indicate which nodes are the duplicated ones
    V(g.sim)$name<-unlist(lapply(V(g.sim)$name,function(x){
        if(x %in% duplicates.names){
            paste("dup-",x,sep="") 
        }else{
            x
        }
    }))
    
    return(g.sim)
}

#########################################################
# Evaluate the clustering results. Take out the duplicates
# and see whether they stay in the same module with the
# original nodes
#########################################################
networkEval<-function(modules){
    
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

get_module<-function(module){
    result<-matrix(ncol=2)
    for(i in 1:nrow(module)){
        module_list<-unlist(strsplit(module[i,2],";"))
        gene<-rep(module[i,1],times=length(module_list))
        new<-cbind(gene,module_list)
        result<-rbind(result,new)
    }
    result<-result[-1,]
    return(result)
}

