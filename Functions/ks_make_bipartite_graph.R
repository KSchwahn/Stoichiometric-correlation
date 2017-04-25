ks_make_bipartite_graph = function(pairs,triplets,quadruples,tr){
library(igraph)
  pairs = pairs[which(pairs$max_correlations>=tr),]
  pairs = pairs[which(pairs$max_pvalues<=0.05),]
  triplets = triplets[which(triplets$adjust_p_value<=0.05),]
  triplets = triplets[which(triplets$correlations>=tr),]
  quadruples = quadruples[which(quadruples$adjust_p_value<=0.05),]
  quadruples = quadruples[which(quadruples$correlations>=tr),]
  nodeSet_A = as.character(pairs$names)
  nodeSet_A_t = paste(sapply(strsplit(as.character(triplets$names),split="[1,2,3,4]\\*"),"[[",2),sapply(strsplit(as.character(triplets$names),split="[1,2,3,4]\\*"),"[[",3),sep="")
  nodeSet_A_q = paste(sapply(strsplit(as.character(quadruples$names),split="[1,2,3,4]\\*"),"[[",2),
                      sapply(strsplit(as.character(quadruples$names),split="[1,2,3,4]\\*"),"[[",3),
                      sapply(strsplit(as.character(quadruples$names),split="[1,2,3,4]\\*"),"[[",4),
                      sapply(strsplit(as.character(quadruples$names),split="[1,2,3,4]\\*"),"[[",5),
                      sep="")
  nodeSet_All <-c(nodeSet_A,nodeSet_A_t,nodeSet_A_q)
  B = lapply(nodeSet_A,function(x){
    tmp_1 = sapply(strsplit(as.character(x),split="->"),"[[",1)
    tmp_2 = sapply(strsplit(as.character(x),split="->"),"[[",2)
    B = c(tmp_1,tmp_2)
  return(B)
  })
  B_t = lapply(nodeSet_A_t, function(x){
    tmp_1 = sapply(strsplit(as.character(x),split="_"),"[[",1)
    tmp_2 = sapply(strsplit(as.character(x),split="_"),"[[",2)
    tmp_2 = sapply(strsplit(as.character(tmp_2),split="->"),"[[",1)
    tmp_3 = sapply(strsplit(as.character(x),split="->"),"[[",2)   
    B_t = c(tmp_1,tmp_2,tmp_3)
    return(B_t)
  })
  B_q = lapply(nodeSet_A_q,function(x){
    tmp_1 = sapply(strsplit(as.character(x),split="_"),"[[",1)
    tmp_2 = sapply(strsplit(as.character(x),split="_"),"[[",2)
    tmp_2 = sapply(strsplit(as.character(tmp_2),split="->"),"[[",1)
    tmp_3 = sapply(strsplit(as.character(x),split="->"),"[[",2)
    tmp_3 = sapply(strsplit(as.character(tmp_3),split="_"),"[[",1)
    tmp_4 = sapply(strsplit(as.character(x),split="_"),"[[",3)
    B_q = c(tmp_1,tmp_2,tmp_3,tmp_4)  
    return(B_q)
  })
B = c(B,B_t,B_q)
nodeSet_Bll = unique(unlist(B))
  A = lapply(nodeSet_A,function(x){
    A = rep(x,2)
    return(A)
  })
  A_t = lapply(nodeSet_A_t,function(x){
	  A_t = rep(x,3)
    return(A_t)
  })
  A_q = lapply(nodeSet_A_q,function(x){
    A_q = rep(x,4)
    return(A_q) 
  })
A = c(A,A_t,A_q)
  g <- graph.empty()
  g <- add.vertices(g,nv=length(nodeSet_All),attr=list(name=as.character(nodeSet_All),type=rep(TRUE,length(nodeSet_All))))
  g <- add.vertices(g,nv=length(nodeSet_Bll),attr=list(name=as.character(nodeSet_Bll),type=rep(FALSE,length(nodeSet_Bll))))

  edgeList = data.frame(A=unlist(A),B=unlist(B))
  # we need to turn edgeList into a vector (and using names instead of indexes)
  edgeListVec <- as.vector(t(as.matrix(data.frame(S1=edgeList$A,S2=edgeList$B))))
  g <- add.edges(g,edgeListVec)
  W = lapply(pairs$max_correlation,function(x){
    W = rep(x,2)
    return(W)
  })
  W_t = lapply(triplets$correlations,function(x){
    W_t = rep(x,3)
    return(W_t)
  })
  W_q = lapply(quadruples$correlations,function(x){
    W_q = rep(x,4)
    return(W_q) 
  })
  Weight<-c(unlist(W),unlist(W_t),unlist(W_q))
  E(g)$weight<-Weight
  print(is.bipartite(g))
return(g)
}


