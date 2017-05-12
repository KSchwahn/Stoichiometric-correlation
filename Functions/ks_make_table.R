ks_make_table = function(pair,Corr,Names,tr){
#The function takes all calculated pairs, triples and quadruples and calculates several statistics, as max, min, mean correlation for each metabolite. It is separated for pairs, triples and quadruples

####Input
  #pairs: supply the output of the functions "ks_find_max_cor()"
  #Corr: supply the result of the "File_read_triples.py-script" and "File_read_quadruples.py-script" in one list. These document has to be read into R first.
  #tr: specify the threshold for which the statistics should be calculated.

####Output
  #OUT: a data.frame that list for each metabolite the following informations
	#-Total_number_of_correlations
	#-Triplet_number_correlation
	#-Quadruple_number_correlation
	#-Pairs_number_correlation
	#-Triplet_mean_correlation
	#-Quadruple_mean_correlation
	#-Pairs_mean_correlation
	#-Triplet_max_correlation
	#-Quadruple_max_correlation
	#-Pairs_max_correlation
	#-Triplet_min_correlation
	#-Quadruple_min_correlation
	#-Pairs_min_correlation
	#-Stoichiometric_mean -> mean of the reported indices 
	#-Stoichiometric_max -> max of the reported indices


  pair = pair[which(pair$max_pvalues<=0.05),]
  trip = Corr$triplets[which(Corr$triplets$adjust_p_value<=0.05),]
  quad = Corr$quadruples[which(Corr$quadruples$adjust_p_value<=0.05),]
  
  pair = pair[which(pair$max_correlations>=tr),]
  trip = trip[which(trip$correlations>=tr),]
  quad = quad[which(quad$correlations>=tr),]
  
  total = c()
  t_num = c()
  q_num = c()
  p_num = c()
  t_mean = c()
  q_mean = c()
  p_mean = c()
  t_max = c()
  q_max = c()
  p_max = c()
  t_min = c()
  q_min = c()
  p_min = c()
  sto_max = c()
  sto_mean = c()
  
  
  for(i in Names){
    print(i)
    ind_t=grep(as.character(i),trip$names)
    ind_q=grep(as.character(i),quad$names)
    ind_p=grep(i,pair$names)
    tmp_trip<-trip[ind_t,]
    tmp_quad<-quad[ind_q,]
    tmp_pair<-pair[ind_p,]
    #
    
    total = c(total,sum(nrow(tmp_trip),nrow(tmp_quad),nrow(tmp_pair)))
    t_num = c(t_num,nrow(tmp_trip))
    q_num = c(q_num,nrow(tmp_quad))
    p_num = c(p_num,nrow(tmp_pair))
    t_mean=c(t_mean,mean(tmp_trip$correlations))
    q_mean=c(q_mean,mean(tmp_quad$correlations))
    p_mean=c(p_mean,mean(tmp_pair$max_correlations))
    t_max=c(t_max,max(tmp_trip$correlations))
    q_max=c(q_max,max(tmp_quad$correlations))
    p_max=c(p_max,max(tmp_pair$max_correlations))
    t_min=c(t_min,min(tmp_trip$correlations))
    q_min=c(q_min,min(tmp_quad$correlations))
    p_min=c(p_min,min(tmp_pair$max_correlations))
    
    if(length(tmp_trip$names)&length(tmp_quad$names)>0){
      
      c_1 = sapply(strsplit(as.character(tmp_trip$names),split="\\*"),"[[",1)
      c_2 = sapply(strsplit(sapply(strsplit(as.character(tmp_trip$names),split="\\*"),"[[",2),split="_"),"[[",2)
      c_3 = sapply(strsplit(as.character(tmp_quad$names),split="\\*"),"[[",1)
      c_4 = sapply(strsplit(sapply(strsplit(as.character(tmp_quad$names),split="\\*"),"[[",2),split="_"),"[[",2)
      c_5 = sapply(strsplit(sapply(strsplit(as.character(tmp_quad$names),split="\\*"),"[[",3),split="->"),"[[",2)
      c_6 = sapply(strsplit(sapply(strsplit(as.character(tmp_quad$names),split="\\*"),"[[",4),split="_"),"[[",2)
      sto_max = c(sto_max,max(c(as.numeric(c_1),as.numeric(c_2),as.numeric(c_3),as.numeric(c_4),as.numeric(c_5),as.numeric(c_6))))
      sto_mean = c(sto_mean,mean(c(as.numeric(c_1),as.numeric(c_2),as.numeric(c_3),as.numeric(c_4),as.numeric(c_5),as.numeric(c_6))))
    }
    else{
      sto_max = c(sto_max,0)
      sto_mean = c(sto_mean,0)
    }
    
  }

  OUT<-data.frame(Names=Names,Total_number_of_correlations=total,Triplet_number_correlation=t_num,Quadruple_number_correlation=q_num,Pairs_number_correlation=p_num,
                  Triplet_mean_correlation=t_mean,Quadruple_mean_correlation=q_mean,Pairs_mean_correlation=p_mean,
                  Triplet_max_correlation=t_max,Quadruple_max_correlation=q_max,Pairs_max_correlation=p_max,
                  Triplet_min_correlation=t_min,Quadruple_min_correlation=q_min,Pairs_min_correlation=p_min,
                  Stoichiometric_mean = sto_mean,Stoichiometric_max = sto_max)
  OUT[OUT==Inf]<-NA
  OUT[OUT==-Inf]<-NA
  return(OUT)
}
