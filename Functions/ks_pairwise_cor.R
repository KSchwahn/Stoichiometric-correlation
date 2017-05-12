####Set of functions to caculate the pairwise correlation of the provided metabolites and to extract the maximal correlation in a way it matches the output of "ks_find_max_cor_tr()" and "ks_find_max_cor_qu()"

ks_pairwise_cor = function(Data,log=T){
####Input
  #Data: data.frame with metabolite measurements. Rows: metabolites; Cols: experiments
  #log: specifies of a logtransformation should be performed as it is done for the triples and quadruples

#####Output
  #List with two data.frames. One with all correlation values, the second with the associated adjusted p-values.

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
###Input
  #Data: The list provided as a output by "ks_pairwise_cor()"

###Output
  #A data.frame, which matches the output of the triples and quadruple analysis. It contains the maximal correlation, the associated p-value and the names of the metabolites, which form the pair.

  correlations = Data$correlations
  pvalues = Data$pvalues
  out_correlations = numeric(length(correlations[upper.tri(correlations,diag=F)]))
  out_pvalues = numeric(length(correlations[upper.tri(correlations,diag=F)]))
  out_names = character(length(correlations[upper.tri(correlations,diag=F)]))
  counter = 0
  for(i in 2:nrow(correlations)-1){

    for(j in (i+1):ncol(correlations)){

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


