ks_find_max_cor_qu = function(quadruples,indices,NRcluster1=1,tr=0.8){
  library(parallel)
  ind = c(combn(1:length(indices),2,simplify=F),combn(length(indices):1,2,simplify=F))
  frame2=length(ind)*length(ind)
 
  t_1=sapply(strsplit(sapply(strsplit(as.character(quadruples$names),split="\\*"),"[[",2),split="_"),"[[",1)
  t_2=sapply(strsplit(sapply(strsplit(as.character(quadruples$names),split="\\*"),"[[",3),split="->"),"[[",1)
  t_3=sapply(strsplit(sapply(strsplit(as.character(quadruples$names),split="\\*"),"[[",4),split="_"),"[[",1)
  t_4=sapply(strsplit(as.character(quadruples$names),split="\\*"),"[[",5)
  c_1 = sapply(strsplit(as.character(quadruples$names),split="\\*"),"[[",1)
  c_2 = sapply(strsplit(sapply(strsplit(as.character(quadruples$names),split="\\*"),"[[",2),split="_"),"[[",2)
  c_3 = sapply(strsplit(sapply(strsplit(as.character(quadruples$names),split="\\*"),"[[",3),split="->"),"[[",2)
  c_4 = sapply(strsplit(sapply(strsplit(as.character(quadruples$names),split="\\*"),"[[",4),split="_"),"[[",2)
  c_5 = as.numeric(c_1)+as.numeric(c_2)+as.numeric(c_3)+as.numeric(c_4)
 
  
 qua_names=paste(paste(t_1,t_2,sep=" "),paste(t_3,t_4,sep=" "),sep="_")

 index = which(quadruples$correlations>=tr & quadruples$adjust_p_value<=0.05)
 quadruples = quadruples[index,]
 qua_names = qua_names[index]
 c_5=c_5[index]

 store = c()
 unique_names = unique(qua_names)
 name_1 = c()
 name_2 = c()
 for(i in 1:length(unique_names)){
   print(i)
   changed = paste(unlist(strsplit(unique_names[i],"_"))[2],unlist(strsplit(unique_names[i],"_"))[1],sep="_")  
   store = c(store,changed)
   if(unique_names[i]%in%store){
     next
   }
   else{
     name_1 = c(name_1,unique_names[i])
     name_2 = c(name_2,changed)
   }
 }
 DF = data.frame(V1 = name_1,V2=name_2)
 rows_DF = nrow(DF)
 RES = mclapply(1:rows_DF,mc.cores=NRcluster1,mc.preschedule=F,function(y){
    print(y)
    tmp = quadruples[grep(paste("^",DF[y,1],sep=""),qua_names),]
    tmp_c = c_5[grep(paste("^",DF[y,1],sep=""),qua_names)]
    
    tmp  = rbind(tmp,quadruples[grep(paste("^",DF[y,2],sep=""),qua_names),])
    tmp_c = c(tmp_c,c_5[grep(paste("^",DF[y,2],sep=""),qua_names)])
    
    #print(dim(tmp))
    #print(DF[y,])
    out=tmp[which(tmp$correlations==max(tmp$correlations)),]
    tmp_c=tmp_c[which(tmp$correlations==max(tmp$correlations))]
    
    if(nrow(out)!=1){
      out=out[which(tmp_c==min(tmp_c)),]
      #print("larger 1")
      write.table(out,file=paste(as.character(out$names),"Q.tab",sep=""),row.names = F,col.names = T,quote = F,sep = "\t")
      
    }
    else{
      #print("is 1")
      write.table(out,file=paste(as.character(out$names),"Q.tab",sep=""),row.names = F,col.names = T,quote = F,sep = "\t")
      
    }
    
  })
  return(NULL)
}

