ks_stoichiometric_correlation = function(Data,indices=c(1:4),nblocks, NRcluster=1,names){
  #Important: Rows of data are metabolites
  #nrows of data %% nblokcs musts gives remainder 0

  #The function calculates all possible combination of triples and quadruples, based on the input data and the indices
  #Up following functions will calculate the maximal correlations for triples and quadruples for each set of metabolites

####Input
   #Data: data.frame with metabolite measurements. Rows: metabolites; Cols: experiments
   #indicies: set of values with which the metabolites are multiplied for the triples and quadruples. The more indices are specified more combinations are calculated. Default: [1,2,3,4]
   #nblocks: Defines the size of the matrix blockes passed to rcorr. The provided function "divisors()" gives the potential nblocks, that can be chosen. As a default, the number of rows of the Data can be given to the function
   #NRcluster: Defines how many cores should be used for the parallelization. The default is 1, which means that the task is not parallelized. The parallelization is performed on the step of correlation calculation and more Blocks can be used simultaneously.
   #names: Separated input of the metabolite names. This is needed for correctly assigning the triples and quadruples and to make the evaluation of the final output simpler.

####Output
   #A large list, which holds all possible combinations of triples and quadruples. The list is likely to be unsorted, based on splitting the data into blocks. There for pass the results to "ks_find_max_cor_tr()" and "ks_find_max_cor_qu()".

  #load needed packages
  library(Hmisc, quietly = TRUE)
  library(igraph, quietly = TRUE)
  library(parallel)
  print("Start calculating the stoichiometric correlations")
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
    print(paste("Reading file",i,sep=" "))
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


