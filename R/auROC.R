#source("clustering.R")
library(igraph)
library(ROCR)
library(ggplot2)

#####################################################
# Compute ROC curve and AUC values
#####################################################

computeAuROC<-function(g,nDup,dDup,rewProb,simMethods,nThreads=1){
    
    # Check type of input argument
    if(!is(g,"igraph")){
        stop("Input graph should be a igraph graph object!")  
    }
    if(!is(nDup,"numeric")){
        stop("nDup should be an numeric value!")
    }
    if(!is(dDup,"numeric")){
        stop("dDup should be an numeric value!")
    }
    if(!is(rewProb,"numeric")){
        stop("rewProb should be an numeric value!")
    }

    # Duplicate the graph  
    g.dup<-networkDup(g,nDup,dDup)
  
    # Rewire the duplicates, get the score for the merged graph
    g.dup.rewired<-networkRew(g.dup,rate=rewProb)
    g.merge<-g+g.dup.rewired
    
    # Compute a ranked list of similarity scores for each similarity measure
    lst.rank<-list()
    for(method in simMethods){
        
        if(method != "invlogweighted" && method != "jaccard" && method != "geometric" && method !="wt"){
            stop("Invalid value for simMethods!")  
        }
      
        lst.rank[[length(lst.rank)+1]]<-switch(    method,
                                                   invlogweighted = get_rank_for_coreg(g.merge,sim="invlogweighted",mode="duplicates"),
                                                   jaccard = get_rank_for_coreg(g.merge,sim="jaccard",mode="duplicates"),
                                                   geometric = get_rank_for_coreg(g.merge,sim="geometric",mode="duplicates",nThreads=nThreads),
                                                   wt = get_rank_for_wt(g.merge)
                                               )
    }
    
    AUC<-c()
    
    # Format the result and feed into functions that caculate the ROC and AUC
    fig<-list()
    curve<-data.frame()
    for(i in 1:length(lst.rank)){
            pred<-prediction(lst.rank[[i]][,3],as.factor(lst.rank[[i]][,4]))
            perf<-performance(pred,measure="tpr",x.measure="fpr")
            auc<-performance(pred,measure="auc")@y.values[[1]]
            AUC<-c(AUC,auc)
            curve<-rbind(curve,data.frame(unlist(perf@x.values),unlist(perf@y.values),rep(simMethods[i],times=length(unlist(perf@x.values)))))
    }
    names(AUC)<-simMethods
    colnames(curve)<-c("x","y","group")
    re<-list(curve=curve,AUC=AUC)
    class(re)<-"CoReg.auROC"
    return(re)
}

print.CoReg.auROC<-function(auROC){
    cat("CoReg.auROC object\n")
}

plot.CoReg.auROC<-function(auROC){
    # Plot ROC curve
    fig<-ggplot()+geom_path(data=auROC$curve,aes(x=x,y=y,color=group))+
    geom_abline(intercept=0,slope=1)+
                   theme(
                             panel.background=element_blank(),
                             panel.border=element_rect(fill=NA,size=1),
                             axis.ticks=element_line(size=1),
                             text=element_text(size=16),
                             axis.text=element_text(size=16),
                             legend.text=element_text(size=16),
                             legend.key=element_blank()
                        )+
                   xlab("False positive rate")+
                   ylab("True positive rate")+
                   guides(color=guide_legend(title="methods"))+
                   scale_color_manual(values=c("green","red","blue","black"))
    print(fig)
    #return(lst.auc)
}

# Get a ranked list for all the pairs
get_rank_for_coreg<-function(g,sim="invlogweighted",mode="duplicates",nThreads){
  
    if(mode!="duplicates" && mode!="all") stop("mode should be either 'duplicates' or 'all'")
    all.names<-V(g)$name

    # Get duplicates gene
    gene.dup<-all.names[grep('dup-',all.names)]
    gene.dup.names<-unlist(strsplit(gene.dup,split="dup-"))
    gene.dup.names<-gene.dup.names[1:length(gene.dup.names)%%2==0]
    
    # Create positive pairs(A->A') and negative pairs('A'->'B')
    gene.pairs.pos<-cbind(gene.dup,gene.dup.names)
    gene.pairs.neg<-t(sapply(1:nrow(gene.pairs.pos),function(x){sample(gene.dup,2,replace=F)},simplify=T))
    gene.pairs<-rbind(gene.pairs.pos,gene.pairs.neg)
    
    # Get similarity matrix. This is only a subset, not a whole network
    if(sim == "invlogweighted"){
        in.sim<-similarity.invlogweighted(g,vids=c(gene.dup,gene.dup.names),mode="in")
        out.sim<-similarity.invlogweighted(g,vids=c(gene.dup,gene.dup.names),mode="out")
        
        # Since vids is ignored by inverselogweighted index, similarity matrix needs to be sliced here
        pos<-sapply(c(gene.dup,gene.dup.names),function(x){which(all.names==x)})
        in.sim<-in.sim[,pos]
        out.sim<-out.sim[,pos]
    }else if(sim == "jaccard"){
        in.sim<-similarity.jaccard(g,vids=c(gene.dup,gene.dup.names),mode="in")
        out.sim<-similarity.jaccard(g,vids=c(gene.dup,gene.dup.names),mode="out")
    }else if(sim == "geometric"){
        in.sim<-getSimMatGeo(g,vids=c(gene.dup,gene.dup.names),mode="in",nThreads=nThreads)
        out.sim<-getSimMatGeo(g,vids=c(gene.dup,gene.dup.names),mode="out",nThreads=nThreads)
    }
    
    sim.mat<-in.sim+out.sim
    colnames(sim.mat)<-c(gene.dup,gene.dup.names)
    rownames(sim.mat)<-c(gene.dup,gene.dup.names)
    
    # Get similarity scores for these pairs and mark their class
    scores<-sim.mat[gene.pairs]
    labels<-c(rep(1,nrow(gene.pairs.pos)),rep(0,nrow(gene.pairs.pos)))
    
    re<-data.frame(gene.pairs,scores,labels)
    colnames(re)<-c("gene1","gene2","similarity_score","label")
    re<-re[order(re[,3],decreasing = T),]
    re
}

get_rank_for_wt<-function(g,mode="duplicates"){
  
    if(mode!="duplicates" && mode!="all") stop("Invalid mode. mode should be 'duplicates' or 'all'")  
  
    el<-as_edgelist(g)
    gene.names<-unique(c(el))
    mat.adj<<-matrix(0,nrow=length(gene.names),ncol=length(gene.names),dimnames=list(gene.names,gene.names)) # mat.adj is accessible to other function
    mat.adj[el]<<-1;diag(mat.adj)<<-1  # Add self loop

    # Generate transition matrix and probability matrix
    mat.ts<-t(apply(mat.adj,1,function(x){d<-sum(x);x<-x/d;return(x)}))
    mat.prob<<-mat.ts%*%mat.ts%*%mat.ts%*%mat.ts  # This matrix is accessible to other function

    # Get distance based on the probability matrix
    if(mode=='duplicates'){
        gene.dup<-gene.names[grep('dup-',gene.names)]
        gene.original<-gene.names[-grep('dup-',gene.names)]
        gene.dup.names<-unlist(strsplit(gene.dup,"dup-"))
        gene.dup.names<-gene.dup.names[1:length(gene.dup.names)%%2==0] # Get the original names back
        gene.pos.pairs<-cbind(gene.dup,gene.dup.names)
        gene.neg.pairs<-t(sapply(1:nrow(gene.pos.pairs),function(x){sample(gene.dup,2,replace=F)},simplify=T))
        gene.names.mat<-rbind(gene.pos.pairs,gene.neg.pairs)
    }else{
        gene.names.mat<-t(combn(gene.names,2))
    }
    
    degree.dist<<-apply(mat.adj,2,sum)                                           # A vector storing the incoming degree for each node
    sim<-1-apply(gene.names.mat,1,function(x){get_distance_for_wt(x[1],x[2])})   # similarity(i,j)=1-distance(i,j)
    
    # Get output
    rank.mat<-data.frame(gene.names.mat,sim,c(rep(1,nrow(gene.pos.pairs)),rep(0,nrow(gene.neg.pairs))),stringsAsFactors = F)
    colnames(rank.mat)<-c("gene1","gene2","similarity_score","label")
    rank.mat<-rank.mat[order(rank.mat[,3],decreasing=T),]
    rm(mat.adj,mat.prob,degree.dist,envir=.GlobalEnv)
    rank.mat
}

get_distance_for_wt<-function(i,j){
    return(sqrt(sum((mat.prob[i,]-mat.prob[j,])^2/degree.dist)))
}
