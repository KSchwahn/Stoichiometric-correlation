ks_find_max_cor_tr = function(triplets,indices=c(1:4),NRcluster1=1,NRcluster2=1){
 #The function is used upfollowing of "ks_stochiometric_correlation()". It will calculate for each set of triples (a*M1 + b*M2 -> M3), the one which satisfied the following conditions:
 #- highest correlation of a*M1 + b*M2 to M3
 #- if more than on combinations of [a,b] gives the same correlation, the one with the smallest sum of [a,b] is choosen. The possible values were chosen in "ks_stochiometric_correlation()" with the indices

######Input
  #triples: output of "ks_stochiometric_correlation()". Use only List$triples.
  #indices: set of values with which the metabolites are multiplied for the triples and quadruples. The more indices are specified more combinations are calculated. Default: [1,2,3,4]
  #NRcluster1 + NRcluster2: Defines how many cores should be used for the parallelization. The default is 1, which means that the task is not parallelized. The parallelization is performed on the step of finding the maximal correlation. Setting NRcluster1=1 and NRcluster2=2 in total 2 core are used. Setting NRcluster1=2 and NRcluster2=2 in total 4 core are used.
  

######Output
  #The function will not give back any variables to the R-workspace. It will creat one file per quadruple. This can than be combined with the provided Python-script.


 library(parallel)
  ind = c(combn(1:length(indices),2,simplify=F),combn(length(indices):1,2,simplify=F))
  frame1=length(ind)
  
  t_1=sapply(strsplit(sapply(strsplit(as.character(triplets$names),split="\\*"),"[[",2),split="_"),"[[",1)
  t_2=sapply(strsplit(sapply(strsplit(as.character(triplets$names),split="\\*"),"[[",3),split="->"),"[[",1)
  t_3=sapply(strsplit(sapply(strsplit(as.character(triplets$names),split="\\*"),"[[",3),split="->"),"[[",2)
  c_1=sapply(strsplit(as.character(triplets$names),split="\\*"),"[[",1)
  c_2=sapply(strsplit(sapply(strsplit(as.character(triplets$names),split="\\*"),"[[",2),split="_"),"[[",2)
  c_3=as.numeric(c_1)+as.numeric(c_2)
  
  trip_names=paste(t_1,t_2,t_3,sep=" ")
  rm(t_1,t_2,t_3,c_1,c_2)
  gc()
  triplets=triplets[order(trip_names),]
  c_3 = c_3[order(trip_names)]
  
  Split_t = split(seq(1:nrow(triplets)), ceiling(seq_along(seq(1:nrow(triplets)))/frame1))
  tmp_split=split(seq(1:length(Split_t)),ceiling(seq_along(seq(1:length(Split_t)))/15))
 Max_triplets=mclapply(tmp_split,mc.cores=NRcluster1,mc.preschedule=F,function(y){
    Sub=Split_t[unlist(y)]
    Max_tmp_triplets = mclapply(Sub,function(x){
      tmp=triplets[x,]
      tmp_c = c_3[x]
      out=tmp[which(tmp$correlations==max(tmp$correlations)),]
      tmp_c=tmp_c[which(tmp$correlations==max(tmp$correlations))]
      if(nrow(out)!=1){

        out=out[which(tmp_c==min(tmp_c)),]
        f.name = gsub(pattern = ">",replacement = "",x=gsub(pattern="\\*",replacement = "",x=out$names))
        write.table(out,file=paste(as.character(f.name),"T.tab",sep=""),row.names = F,col.names = T,quote = F,sep = "\t")

      }
      else{
        f.name = gsub(pattern = ">",replacement = "",x=gsub(pattern="\\*",replacement = "",x=out$names))
        write.table(out,file=paste(as.character(f.name),"T.tab",sep=""),row.names = F,col.names = T,quote = F,sep = "\t")
       
      }  
    },mc.cores=NRcluster2,mc.preschedule=F)
    
  gc()
 })
  return(NULL)
}

