ks_stoichiometric_correlation = function(Data,indices,nblocks, NRcluster=1,names){
  #Important: Rows of data are metabolites
  #nrows of data %% nblokcs musst gives remainder 0
  #load needed packages
  library(Hmisc, quietly = TRUE)
  library(igraph, quietly = TRUE)
  library(parallel)
  
  #log transform all data
  Data=log(Data)
  #creat all combinations of length 2 of the metabolite names
  List <- combn(rownames(Data),2,simplify=F)
  #creat all combinations of lenght 2 of the used indices
  ind = c(combn(1:length(indices),2,simplify=F),combn(length(indices):1,2,simplify=F))
  l_mat = c()
  #for each of the combination of data and indices make the sum of the data vectors and creat the matching names
  for(i in 1:length(ind)){
    tmp=t(sapply(List,function(L,i){as.numeric(ind[[i]][1])*Data[L[1],]+as.numeric(ind[[i]][2])*Data[L[2],]},i))
    rownames(tmp)<-unlist(lapply(List,function(L){paste(paste(ind[[i]][1],L[1],sep="*"),paste(ind[[i]][2],L[2],sep="*"),sep="_")}))
    l_mat = rbind(l_mat,tmp)
  }
  #pass the complete matrix over to the modified bigcor function
  #the result is a large data frame with all correlations and p-values
  #Important: nblocks defines the size of the matrix blockes passed to rcorr
  #it musst be an integer so that ncol(t(L_mat)) %% nblocks == 0
  L_mat = rbind(Data,l_mat)
  print(nrow(L_mat))
  G= bigcor_p(t(L_mat),nblocks=nblocks,type="pearson",names, NRcluster)
  G_list = c()
  files=list.files("./","Block*")
  for(i in 1:length(files)){
    print(i)
    G<-read.table(files[i],header = T,sep = "\t",stringsAsFactors = F)
    G_list = rbind(G_list,G)
    rm(G)
    gc()
  }
  G <- G_list
  rm(G_list)
  gc()
  names=paste(G$to,G$from,sep="->")
  G = cbind(G,names)
  #split the data frame into triplets and quadruples, pair correlaton can be easily computed with the normal rcorr function
  G_q<-G[-which(G$from%in%rownames(Data)|G$to%in%rownames(Data)),]
  G_t<-G[!(G$names%in%G_q$names),]
  
  G_t = data.frame(names=as.character(G_t$names),correlations=as.numeric(G_t$weight),p_value=as.numeric(G_t$p_value),
                   adjust_p_value=as.numeric(p.adjust(G_t$p_value,method = "BH")))
  G_t = unique(G_t)
  
  
  G_q = data.frame(names=as.character(G_q$names),correlations=as.numeric(G_q$weight),p_value=as.numeric(G_q$p_value),
                   adjust_p_value=as.numeric(p.adjust(G_q$p_value,method = "BH")))
  
  G_q = unique(G_q)
  #return the results in a list
  return(list(triplets=G_t,quadruples = G_q))
}


