ks_shared_metabolites = function(data_pair,data_all,data2_pair,data2_all,tr,name1,name2){
#This functions is important when comparing two data sets.
#It computes the overlap of pairs, triples and quadruples of two data sets

###Input
 #data_pair: supply the output of the functions "ks_find_max_cor()" for the first data set
 #data_pair2: supply the output of the functions "ks_find_max_cor()" for the second data set
 #data_all: supply the result of the "File_read_triples.py-script" and "File_read_quadruples.py-script" for the first data set in one list. These document has to be read into R first.
 #data_all2: supply the result of the "File_read_triples.py-script" and "File_read_quadruples.py-script" for the second data set in one list. These document has to be read into R first.
 #tr: specify the threshold for which the soverlap should be calculated.
 #name1: name of the first data set. Supply a string, which will be used for renaming the output.
 #name2: name of the second data set. Supply a string, which will be used for renaming the output.

###Output
 #The function returns a list with four data.frame objects.
 #Counts: overview over the overlapping pairs, triples and quadruples and the overall number per data set.
 #Pairs: On column with the names of the overlapping pairs.
 #Triples: On column with the names of the overlapping triples.
 #Quadruples: On column with the names of the overlapping quadruples.


  data_pair = data_pair[which(data_pair$max_pvalues<=0.05 & data_pair$max_correlations>=tr),]
  data_trip = data_all$triplets[which(data_all$triplets$adjust_p_value<=0.05 & data_all$triplets$correlations>=tr),]
  data_quad = data_all$quadruples[which(data_all$quadruples$adjust_p_value<=0.05 & data_all$quadruples$correlations>=tr),]
  
  data2_pair = data2_pair[which(data2_pair$max_pvalues<=0.05 & data2_pair$max_correlations>=tr),]
  data2_trip = data2_all$triplets[which(data2_all$triplets$adjust_p_value<=0.05 & data2_all$triplets$correlations>=tr),]
  data2_quad = data2_all$quadruples[which(data2_all$quadruples$adjust_p_value<=0.05 & data2_all$quadruples$correlations>=tr),]
  
  
  data_t_names=paste(sapply(strsplit(as.character(data_trip$names),split="[1,2,3,4]\\*"),"[[",2),sapply(strsplit(as.character(data_trip$names),split="[1,2,3,4]\\*"),"[[",3),sep="")
  data_q_names=paste(sapply(strsplit(as.character(data_quad$names),split="[1,2,3,4]\\*"),"[[",2),
                      sapply(strsplit(as.character(data_quad$names),split="[1,2,3,4]\\*"),"[[",3),
                      sapply(strsplit(as.character(data_quad$names),split="[1,2,3,4]\\*"),"[[",4),
                      sapply(strsplit(as.character(data_quad$names),split="[1,2,3,4]\\*"),"[[",5),
                      sep="")
  
  data2_t_names=paste(sapply(strsplit(as.character(data2_trip$names),split="[1,2,3,4]\\*"),"[[",2),sapply(strsplit(as.character(data2_trip$names),split="[1,2,3,4]\\*"),"[[",3),sep="")
  data2_q_names=paste(sapply(strsplit(as.character(data2_quad$names),split="[1,2,3,4]\\*"),"[[",2),
                            sapply(strsplit(as.character(data2_quad$names),split="[1,2,3,4]\\*"),"[[",3),
                            sapply(strsplit(as.character(data2_quad$names),split="[1,2,3,4]\\*"),"[[",4),
                            sapply(strsplit(as.character(data2_quad$names),split="[1,2,3,4]\\*"),"[[",5),
                            sep="") 
  
  
  data_p_names = lapply(strsplit(as.character(data_pair$names),split="->"),sort)
  data2_p_names = lapply(strsplit(as.character(data2_pair$names),split="->"),sort)
  
  Count_pair_data_data2 = length(intersect(data_p_names,data2_p_names))
  
  
  Intersect_pair_data_data2 = intersect(data_p_names,data2_p_names)
  Intersect_pair_data_data2 = unlist(lapply(Intersect_pair_data_data2,function(x){paste(unlist(x),collapse="->")}))
  if(Count_pair_data_data2==0){Intersect_pair_data_data2 = "0"}
    
  data_t_names_2= lapply(sapply(strsplit(data_t_names,split="->"),"[[",2),sort)
  data_t_names_1= lapply(strsplit(sapply(strsplit(data_t_names,split="->"),"[[",1),split="_"),sort)
  data_t_names = paste(data_t_names_1,data_t_names_2)
  data_t_names = mapply(c,data_t_names_1,data_t_names_2, SIMPLIFY=F)
  data_t_names = lapply(data_t_names,function(x){paste(unlist(x),collapse="_")})
  
  data2_t_names_2= lapply(sapply(strsplit(data2_t_names,split="->"),"[[",2),sort)
  data2_t_names_1= lapply(strsplit(sapply(strsplit(data2_t_names,split="->"),"[[",1),split="_"),sort)
  data2_t_names = mapply(c,data2_t_names_1,data2_t_names_2, SIMPLIFY=F)
  data2_t_names = lapply(data2_t_names,function(x){paste(unlist(x),collapse="_")})
  
  Count_triplets_data_data2 = length(intersect(data_t_names,data2_t_names))
  
  
  Intersect_triplets_data_data2 = unlist(intersect(data_t_names,data2_t_names))
  Intersect_triplets_data_data2 = sub("_","+",Intersect_triplets_data_data2)
  Intersect_triplets_data_data2 = sub("_","->",Intersect_triplets_data_data2)
  if(Count_triplets_data_data2==0){Intersect_triplets_data_data2 = "0"}
  
  data_q_names_2= lapply(strsplit(sapply(strsplit(data_q_names,split="->"),"[[",2),split="_"),sort)
  data_q_names_1= lapply(strsplit(sapply(strsplit(data_q_names,split="->"),"[[",1),split="_"),sort)
  data_q_names = mapply(c,data_q_names_1,data_q_names_2, SIMPLIFY=F)
  data_q_names = lapply(data_q_names,function(x){paste(unlist(x),collapse="_")})
  
  data2_q_names_2= lapply(strsplit(sapply(strsplit(data2_q_names,split="->"),"[[",2),split="_"),sort)
  data2_q_names_1= lapply(strsplit(sapply(strsplit(data2_q_names,split="->"),"[[",1),split="_"),sort)
  data2_q_names = mapply(c,data2_q_names_1,data2_q_names_2, SIMPLIFY=F)
  data2_q_names = lapply(data2_q_names,function(x){paste(unlist(x),collapse="_")})
  
  Count_quadruples_data_data2 = length(intersect(data_q_names,data2_q_names))
  
  Intersect_quadruples_data_data2 = unlist(intersect(data_q_names,data2_q_names))
  Intersect_quadruples_data_data2 = sub("_","+",Intersect_quadruples_data_data2)
  Intersect_quadruples_data_data2 = sub("_","->",Intersect_quadruples_data_data2)
  Intersect_quadruples_data_data2 = sub("_","+",Intersect_quadruples_data_data2)
  if(Count_quadruples_data_data2==0){Intersect_quadruples_data_data2 = "0"}
  
  #rename the output
  str1 = paste("Shared between",name1,"and",name2,sep=" ")
  str2 = name1
  str3 = name2
  
  return(list(Counts=data.frame(Compare=c(str2,str3,str1),
                                Pairs=c(length(data_p_names),
                                        length(data2_p_names),
                                        Count_pair_data_data2),
                                Triplets=c(length(data_t_names),
                                           length(data2_t_names),
                                           Count_triplets_data_data2),
                                Quadruples=c(length(data_q_names),
                                             length(data2_q_names),
                                             Count_quadruples_data_data2)),
              pair_data_data2 = Intersect_pair_data_data2,
              triplets_data_data2 = Intersect_triplets_data_data2,
              quadruples_data_data2 = Intersect_quadruples_data_data2))
  
}
