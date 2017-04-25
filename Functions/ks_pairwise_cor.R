ks_pairwise_cor = function(Data,log=T){
library(Hmisc)
  if(log==T){
    corr = rcorr(log(t(Data)),type="pearson")
    mat<-corr$r
    p_value<-corr$P
    p_value=matrix(p.adjust(p_value,method="BH"),nrow=nrow(Data),ncol=nrow(Data),byrow=T)
   
  }
  else{
    corr = rcorr(t(Data),type="pearson")
    mat<-corr$r
    p_value<-corr$P
    p_value=matrix(p.adjust(p_value,method="BH"),nrow=nrow(Data),ncol=nrow(Data),byrow=T)
   
  }
  rownames(mat)<-rownames(Data)
  colnames(mat)<-rownames(Data)
  rownames(p_value)<-rownames(Data)
  colnames(p_value)<-rownames(Data)
  diag(mat)<-0
   
  return(list(correlations = mat, pvalues=p_value))
}


ks_find_max_cor = function(Data){
  correlations = Data$correlations
  pvalues = Data$pvalues
  out_correlations = numeric(length(correlations[upper.tri(correlations,diag=F)]))
  out_pvalues = numeric(length(correlations[upper.tri(correlations,diag=F)]))
  out_names = character(length(correlations[upper.tri(correlations,diag=F)]))
  counter = 0
  for(i in 2:nrow(correlations)-1){
    print(i)
    for(j in (i+1):ncol(correlations)){
      print(j)
      if(i!=j){
        counter = counter+ 1 
        out_correlations[counter]<-correlations[i,j]
        out_pvalues[counter]<-pvalues[i,j]
        out_names[counter]<-paste(rownames(correlations)[i],colnames(correlations)[j],sep="->")
      }
    }
  }
  
  return(data.frame(max_correlations=out_correlations,max_pvalues=out_pvalues,names=out_names))
}


