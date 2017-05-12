bigcor_p <- function(x, nblocks = 10, type = "pearson", names, NRcluster){
  #modified function from http://www.r-bloggers.com/bigcor-large-correlation-matrices-in-r/
  #This function is called by "ks_stoichiometric_correlation". 
  #It will calculate the correlation matrices and write the files containing these into the working directory. 
  NCOL <- ncol(x)
  ## test if ncol(x) %% nblocks gives remainder 0
  if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
  ## split column numbers into 'nblocks' groups
  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  ## create all unique combinations of blocks
  COMBS <- combn(1:length(SPLIT), 2,simplify = F)
  Glist <- mclapply(COMBS, function(COMB){
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    print(paste("Block", COMB[1], "with Block", COMB[2], sep=" "))
    flush.console()
    #use rcorr instead of cor, for the p-values and speed gain
    COR <- rcorr(x[, G1], x[, G2], type=type)  

    #remove na and add value to the p-values so it is not lost during the weighting of the igraph package
    COR$P[is.na(COR$P)] <- 0
    COR$P <- COR$P + 0.001
    diag(COR$r) <- 0
    diag(COR$P) <- 0
    if(any(names%in%rownames(COR$r))){
      COR$r[names,names] <- 0
      COR$P[names,names] <- 0
    }
    G_w = get.data.frame(graph.adjacency(COR$r,diag=F,weighted = T,mode = "undirected"))
    G_p = get.data.frame(graph.adjacency(COR$P,diag=T,weighted = T,mode = "undirected"))
    G_tmp <- cbind(G_w,p_value=G_p$weight-0.001) 
    G = G_tmp
    rm(G_w,G_p,G_tmp)
    G = G[!duplicated(G),]
    gc()
    G=G[!apply(G,1,function(R){a<-unlist(sapply(unlist(strsplit(R[1],split="_")),function(x){strsplit(x,split="[0-9]\\*")}));
                               b<-unlist(sapply(unlist(strsplit(R[2],split="_")),function(x){strsplit(x,split="[0-9]\\*")}));
                               any(a[which(a!="")]%in%b[which(b!="")])}),]
    
    COR <- NULL
    gc()
    print("Writing to file")
    write.table(G,file=paste("Block", COMB[1], "with Block", COMB[2],".tab",sep=""),quote = F,row.names = F,col.names = T,sep="\t")
    rm(G)
    gc()
  },mc.cores=NRcluster)
  
  
  gc()
  return(NULL)
}
