
# The example bla bla bla... please explain what is this example and what you get as a result!

# Preparation Step

# loading data:
# metabolites in rows and samples/conditions in columns

# reading A. thaliana's metabolite data set
Ara = read.table("Ara_data.txt",sep="\t",header = T)
# reading E. Coli's metabolite data set
Ecoli = read.table("Ecoli_data.txt",sep="\t",header = T)

# sourcing R-functions in the corresponding folder folder
file.sources = list.files(path = "../Functions",pattern="/*.R",full.names = T)
sapply(file.sources,source,.GlobalEnv)

# creating temporary folders before running SCA
subDir = "Ecoli_run/"

if (file.exists(subDir)){
  setwd(subDir)
} else {
  dir.create(subDir)
  setwd(subDir)
}

# Setting input arguments before calling the function ks_stoichiometric_correlation which calculates the stoichiometric correlations (SCs).
# Data is the metabolomics data set.
# indices is assigned to the set of {1,2,3,4} as explained in the manuscript. 
# nblocks, Number of blocks, is set to the number of metabolites. This number will be passed to the function divirors which identifies the 
# number of blocks into which the data can be divided. Please refer to the user's manual or manuscript for more details. 
# NRcluster shows the number of available nodes for calculating the SCs.
# Set NRcluster to a value larger than 1 if you have the possibility to run the program in parallel mode.
NRcluster <- 1
# names, the name of metabolites should be passed to the function

# Calculates the SCs for all pairs, triplets, and quadruples. 
Ecoli_cor = ks_stoichiometric_correlation(Data = Ecoli,indices = c(1:4), nblocks = nrow(Ecoli), NRcluster = NRcluster, names = rownames(Ecoli))


Ecoli_tr=ks_find_max_cor_tr(triplets = Ecoli_cor$triplets, indices = c(1:4),NRcluster1 = 1, NRcluster2 = 1)

Ecoli_qu=ks_find_max_cor_qu(Ecoli_cor$quadruples,c(1:4),NRcluster,tr=0.8)

#run python scripts
command = "python"
path2trip='"../../Functions/File_read_triples.py"'
path2quad='"../../Functions/File_read_quadruples.py"'

# Add path to script as first arg
allArgs = c(path2trip,'"Ecoli_triplets.txt"')
output = system2(command, args=allArgs, stdout=TRUE)
allArgs = c(path2quad,'"Ecoli_quadruples.txt"')
output = system2(command, args=allArgs, stdout=TRUE)

Ecoli_maxcor_tr = read.table("Ecoli_triplets.txt",header=T,sep = "\t")
Ecoli_maxcor_qu = read.table("Ecoli_quadruples.txt",header=T,sep="\t")

Ecoli_MaxCor = list(triplets=Ecoli_maxcor_tr,quadruples=Ecoli_maxcor_qu)

#Calculate pairwise correlation between metabolites. This is similar to standard pearson correlation.
#Call log=F if no log transformation is wanted
Ecoli_pair_log=ks_pairwise_cor(Ecoli,log=T)
Ecoli_pair_log_max=ks_find_max_cor(Ecoli_pair_log)

Ecoli_pair=ks_pairwise_cor(Ecoli,log=F)
Ecoli_pair_max=ks_find_max_cor(Ecoli_pair)

setwd("../")
subDir = "Ara_run/"

if (file.exists(subDir)){
  setwd(subDir)
} else {
  dir.create(subDir)
  setwd(subDir)
}

#increase NRcluster if more cores are available on your system. NRcluster=1 means that the task is not running in parallel.
Ara_cor = ks_stoichiometric_correlation(Ara,c(1:4),nrow(Ara), NRcluster= NRcluster,rownames(Ara))

Ara_tr=ks_find_max_cor_tr(triplets=Ara_cor$triplets,c(1:4),1,NRcluster)

Ara_qu=ks_find_max_cor_qu(Ara_cor$quadruples,c(1:4),NRcluster,tr=0.8)

allArgs = c(path2trip,'"Ara_triplets.txt"')
output = system2(command, args=allArgs, stdout=TRUE)
allArgs = c(path2quad,'"Ara_quadruples.txt"')
output = system2(command, args=allArgs, stdout=TRUE)

#Load than the files and remove all intermediate files if needed.
Ara_maxcor_tr = read.table("Ara_triplets.txt",header=T,sep = "\t")
Ara_maxcor_qu = read.table("Ara_quadruples.txt",header=T,sep="\t")

Ara_MaxCor = list(triplets=Ara_maxcor_tr,quadruples=Ara_maxcor_qu)

#Calculate pairwise correlation between metabolites. This is similar to standard pearson correlation.
#Call log=F if no log transformation is wanted
Ara_pair_log=ks_pairwise_cor(Ara,log=T)
Ara_pair_log_max=ks_find_max_cor(Ara_pair_log)

Ara_pair=ks_pairwise_cor(Ara,log=F)
Ara_pair_max=ks_find_max_cor(Ara_pair)


setwd("../")

##########################################################
write.table(ks_make_table(Ecoli_pair_log_max,Ecoli_MaxCor,rownames(Ecoli),0.8),file = "Ecoli_complete_table_08.tab",row.names=F,col.names=T,quote = F,sep="\t")
write.table(ks_make_table(Ecoli_pair_log_max,Ecoli_MaxCor,rownames(Ecoli),0.85),file = "Ecoli_complete_table_085.tab",row.names=F,col.names=T,quote = F,sep="\t")
write.table(ks_make_table(Ecoli_pair_log_max,Ecoli_MaxCor,rownames(Ecoli),0.9),file = "Ecoli_complete_table_09.tab",row.names=F,col.names=T,quote = F,sep="\t")

write.table(ks_make_table(Ara_pair_log_max,Ara_MaxCor,rownames(Ara),0.8),file = "Ara_complete_table_08.tab",row.names=F,col.names=T,quote = F,sep="\t")
write.table(ks_make_table(Ara_pair_log_max,Ara_MaxCor,rownames(Ara),0.85),file = "Ara_complete_table_085.tab",row.names=F,col.names=T,quote = F,sep="\t")
write.table(ks_make_table(Ara_pair_log_max,Ara_MaxCor,rownames(Ara),0.9),file = "Ara_complete_table_09.tab",row.names=F,col.names=T,quote = F,sep="\t")


############################################################
#Procede from here if you want to compare two data sets:

Ecoli_bipartite_graph_08 = ks_make_bipartite_graph(Ecoli_pair_log_max,Ecoli_MaxCor$triplets,Ecoli_MaxCor$quadruples,0.8)
Ecoli_bipartite_graph_085 = ks_make_bipartite_graph(Ecoli_pair_log_max,Ecoli_MaxCor$triplets,Ecoli_MaxCor$quadruples,0.85)
Ecoli_bipartite_graph_09 = ks_make_bipartite_graph(Ecoli_pair_log_max,Ecoli_MaxCor$triplets,Ecoli_MaxCor$quadruples,0.9)

Ara_bipartite_graph_08 = ks_make_bipartite_graph(Ara_pair_log_max,Ara_MaxCor$triplets,Ara_MaxCor$quadruples,0.8)
Ara_bipartite_graph_085 = ks_make_bipartite_graph(Ara_pair_log_max,Ara_MaxCor$triplets,Ara_MaxCor$quadruples,0.85)
Ara_bipartite_graph_09 = ks_make_bipartite_graph(Ara_pair_log_max,Ara_MaxCor$triplets,Ara_MaxCor$quadruples,0.9)


data_08= data.frame(names(sort(degree(Ecoli_bipartite_graph_08,v=rownames(Ecoli)))),sort(degree(Ecoli_bipartite_graph_08,v=rownames(Ecoli))),
                      names(sort(degree(Ara_bipartite_graph_08,v=rownames(Ara)))),sort(degree(Ara_bipartite_graph_08,v=rownames(Ara))),
                      row.names = NULL)
colnames(data_08)<-c("Metabolites","Ecoli metabolite coupling degree","Metabolites","Ara metabolite coupling degree")

data_085= data.frame(names(sort(degree(Ecoli_bipartite_graph_085,v=rownames(Ecoli)))),sort(degree(Ecoli_bipartite_graph_085,v=rownames(Ecoli))),
                    names(sort(degree(Ara_bipartite_graph_085,v=rownames(Ara)))),sort(degree(Ara_bipartite_graph_085,v=rownames(Ara))),
                    row.names = NULL)
colnames(data_085)<-c("Metabolites","Ecoli metabolite coupling degree","Metabolites","Ara metabolite coupling degree")

data_09 = data.frame(names(sort(degree(Ecoli_bipartite_graph_09,v=rownames(Ecoli)))),sort(degree(Ecoli_bipartite_graph_09,v=rownames(Ecoli))),
                    names(sort(degree(Ara_bipartite_graph_09,v=rownames(Ara)))),sort(degree(Ara_bipartite_graph_09,v=rownames(Ara))),
                    row.names = NULL)
colnames(data_09)<-c("Metabolites","Ecoli metabolite coupling degree","Metabolites","Ara metabolite coupling degree")

out_list = list(data_08,data_085,data_09)

names(out_list) <- c("Comparison_08","Comparison_085","Comparison_09")
write_list(out_list,wb_name = "Metabolite coupling_degree.xlsx")

###############################

Shared_metabolites_08 = ks_shared_metabolites(Ecoli_pair_log_max,Ecoli_MaxCor,Ara_pair_log_max,Ara_MaxCor,0.8,"Ecoli","Ara")
Shared_metabolites_085 = ks_shared_metabolites(Ecoli_pair_log_max,Ecoli_MaxCor,Ara_pair_log_max,Ara_MaxCor,0.85,"Ecoli","Ara")
Shared_metabolites_09 = ks_shared_metabolites(Ecoli_pair_log_max,Ecoli_MaxCor,Ara_pair_log_max,Ara_MaxCor,0.9,"Ecoli","Ara")

write_list(Shared_metabolites_08,wb_name ="Shared_metabolites_08.xls")
write_list(Shared_metabolites_085,wb_name ="Shared_metabolites_085.xls")
write_list(Shared_metabolites_09,wb_name ="Shared_metabolites_09.xls")

#uncomment to save analysis
save.image("Run.RData")


