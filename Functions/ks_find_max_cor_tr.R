ks_find_max_cor_tr = function(triplets,indices,NRcluster1=1,NRcluster2=1){
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

        write.table(out,file=paste(as.character(out$names),"T.tab",sep=""),row.names = F,col.names = T,quote = F,sep = "\t")

      }
      else{
        write.table(out,file=paste(as.character(out$names),"T.tab",sep=""),row.names = F,col.names = T,quote = F,sep = "\t")
       
      }  
    },mc.cores=NRcluster2,mc.preschedule=F)
    
  gc()
 })
  return(NULL)
}

