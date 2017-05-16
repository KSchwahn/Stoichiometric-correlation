
# This example reproduces the compraison of A.thalina and E. coli from the manuscript 
# "Stoichiometric correlation analysis: principles of metabolic functionality from metabolomics data."

# It loads the meatabolomics data set of both species, calculates the stoichiometric correlations, 
# and finally the the maximal correlations for all possible pairs, triplets and quadruples
# During the calculation temporary folders and files will be created
# Python files are provided to summarize the temporary files and to create files which are loaded back into the R workspace

# Finally, the following output files are generated
# A table including the number of stoichiometric correlations due to pairs, triplets, and quadruples 
# as well as the mean and max of the stoichiometric correlation for each metabolite
# A table comparing the coupling degree of metabolites of both species
# A list of common maximal stoichiometric correlations due to pairs, triplets, and quadruples between the two species

# Preparation Step

# loading data:
# metabolites in rows and conditions/time points in columns

# reading A. thaliana's metabolite data set
Ara = read.table("Ara_data.txt",sep="\t",header = T)
# reading E. Coli's metabolite data set
Ecoli = read.table("Ecoli_data.txt",sep="\t",header = T)

# sourcing R-functions
file.sources = list.files(path = "../Functions",pattern="/*.R",full.names = T)
sapply(file.sources,source,.GlobalEnv)

# creating temporary folder before running SCA
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
# nblocks, Number of blocks, is set to the number of metabolites.
# NRcluster shows the number of available nodes for calculating the SCs.
# Set NRcluster to a value larger than 1 if you have the possibility to run the program in parallel mode.
NRcluster <- 1
# names, the name of metabolites should be passed to the function

# Calculates the SCs for all triplets, and quadruples. 
Ecoli_cor = ks_stoichiometric_correlation(Data = Ecoli,indices = c(1:4), nblocks = nrow(Ecoli), NRcluster = NRcluster, names = rownames(Ecoli))

# Setting input arguments before calling the function ks_find_max_cor_tr which estimates all maximal correlation for triplets.
# triplets are the triplets calculated by ks_stoichiometric_correlation
# indices is assigned to the set of {1,2,3,4} as explained in the manuscript. 
# NRcluster1 and NRcluster2 show the number of available nodes for calculating maximal SCs. 
# The total nodes used for this analysis is the product of NRcluster1 and NRcluster2
# The function will creat temporary files instead of an R-workspace variable
ks_find_max_cor_tr(triplets = Ecoli_cor$triplets, indices = c(1:4),NRcluster1 = 1, NRcluster2 = 1)


# Setting input arguments before calling the function ks_find_max_cor_qu which estimates all maximal correlation for quadruples.
# quadruples are the quadruples calculated by ks_stoichiometric_correlation
# indices is assigned to the set of {1,2,3,4} as explained in the manuscript. 
# NRcluster shows the number of available nodes for calculating the SCs
# tr is the correlation threshold, the quadruples which satisfy the threshold are considered for detecting the maximal correlations. 
# A value above 0 decreases the running time
# The function will creat temporary files instead of an R-workspace variable
ks_find_max_cor_qu(quadruples=Ecoli_cor$quadruples,indices = c(1:4),NRcluster1 = NRcluster,tr=0.8)

# Preparing to load and combine all temporary files with the provided Python scripts
command = "python"

# Define the path to the triplet Python script
path2trip='"../../Functions/File_read_triples.py"'
# Combine the path and name the file to be created by the Python script
allArgs = c(path2trip,'"Ecoli_triplets.txt"')
# Run the Python script
# command calls Python within the terminal
# args is the combined path and file name
output = system2(command, args=allArgs, stdout=TRUE)

# Read the created file into a R-workspace variable
Ecoli_maxcor_tr = read.table("Ecoli_triplets.txt",header=T,sep = "\t")

##
# Define the path to the quadruple Python script
path2quad='"../../Functions/File_read_quadruples.py"'
# Combine the path and name the file to be created by the Python script
allArgs = c(path2quad,'"Ecoli_quadruples.txt"')
# Run the Python script
# command calls Python within the terminal
# args is the combined path and file name
output = system2(command, args=allArgs, stdout=TRUE)

# Read the created file into an R-workspace object
Ecoli_maxcor_qu = read.table("Ecoli_quadruples.txt",header=T,sep="\t")

# Combine the estimated maximal correlations due to triplets and quadruples in one list. 
Ecoli_MaxCor = list(triplets=Ecoli_maxcor_tr,quadruples=Ecoli_maxcor_qu)

#####
# Calculate pairwise correlation between metabolites. 
# Data is the metabolomics data set.
# log has to be set to TRUE to calculate the stoichiometric correlation
Ecoli_pair_log=ks_pairwise_cor(Data = Ecoli,log=T)

# Reformate the output for following analysis
# Data is a list provided as an output of the function ks_pairwise_cor()
Ecoli_pair_log_max=ks_find_max_cor(Data = Ecoli_pair_log)

# Calcualte standard Pearson correlation.
# Data is the metabolomics data set.
Ecoli_pair=ks_pairwise_cor(Data = Ecoli,log=F)

# Reformate the output for upfollowing analysis
# Data is a list provided as an output of the function ks_pairwise_cor()
Ecoli_pair_max=ks_find_max_cor(Data = Ecoli_pair)

##########################################################
#repeat the above described steps for the second data set

setwd("../")
# creating temporary folders before running SCA
subDir = "Ara_run/"

if (file.exists(subDir)){
  setwd(subDir)
} else {
  dir.create(subDir)
  setwd(subDir)
}

# Setting input arguments before calling the function ks_stoichiometric_correlation which calculates the stoichiometric correlations (SCs).
# Data is the metabolomics data set.
# indices is assigned to the set of {1,2,3,4} as explained in the manuscript. 
# nblocks, Number of blocks, is set to the number of metabolites.
# NRcluster shows the number of available nodes for calculating the SCs.
# Set NRcluster to a value larger than 1 if you have the possibility to run the program in parallel mode.
NRcluster <- 1
# names, the name of metabolites should be passed to the function

# Calculates the SCs for all triplets, and quadruples. 
Ara_cor = ks_stoichiometric_correlation(Data=Ara,indices=c(1:4),nblocks=nrow(Ara), NRcluster= NRcluster,names=rownames(Ara))

# Setting input arguments before calling the function ks_find_max_cor_tr which estimates all maximal correlation for triplets.
# triplets are the triplets calculated by ks_stoichiometric_correlation
# indices is assigned to the set of {1,2,3,4} as explained in the manuscript. 
# NRcluster1 and NRcluster2 show the number of available nodes for calculating maximal SCs. 
# The total nodes used for this analysis is the product of NRcluster1 and NRcluster2
# The function will creat temporary files instead of an R-workspace variable
ks_find_max_cor_tr(triplets=Ara_cor$triplets,indices=c(1:4),NRcluster1=1,NRcluster2 = 1)

# Setting input arguments before calling the function ks_find_max_cor_qu which estimates all maximal correlation for quadruples.
# quadruples are the quadruples calculated by ks_stoichiometric_correlation
# indices is assigned to the set of {1,2,3,4} as explained in the manuscript. 
# NRcluster shows the number of available nodes for calculating the SCs
# tr is the correlation threshold, the quadruples which satisfy the threshold are considered for detecting the maximal correlations. 
# A value above 0 decreases the running time
# The function will creat temporary files instead of an R-workspace variable
ks_find_max_cor_qu(quadruples=Ara_cor$quadruples,indices=c(1:4),NRcluster,tr=0.8)


# Preparing to load and combine all temporary files with the provided Python scripts
command = "python"

# Define the path to the triplet Python script
path2trip='"../../Functions/File_read_triples.py"'
# Combine the path and name the file to be created by the Python script
allArgs = c(path2trip,'"Ara_triplets.txt"')

# Run the Python script
# command calls Python within the terminal
# args is the combined path and file name
output = system2(command, args=allArgs, stdout=TRUE)

# Read the created file into an R-workspace object
Ara_maxcor_tr = read.table("Ara_triplets.txt",header=T,sep = "\t")

# Define the path to the quadruple Python script
path2quad='"../../Functions/File_read_quadruples.py"'
# Combine the path and name the file to be created by the Python script
allArgs = c(path2quad,'"Ara_quadruples.txt"')

# Run the Python script
# command calls Python within the terminal
# args is the combined path and file name
output = system2(command, args=allArgs, stdout=TRUE)
# Read the created file into an R-workspace variable
Ara_maxcor_qu = read.table("Ara_quadruples.txt",header=T,sep="\t")

# Combine the estimated maximal correlations due to triplets and quadruples in one list. 
Ara_MaxCor = list(triplets=Ara_maxcor_tr,quadruples=Ara_maxcor_qu)

#####
# Calculate pairwise correlation between metabolites. 
# Data is the metabolomics data set.
# log has to be set to TRUE to calculate the stoichiometric correlation
Ara_pair_log=ks_pairwise_cor(Data=Ara,log=T)

# Reformate the output for following analysis
# Data is a list provided as an output of the function ks_pairwise_cor()
Ara_pair_log_max=ks_find_max_cor(Data=Ara_pair_log)

# Calcualte standard Pearson correlation.
# Data is the metabolomics data set.
Ara_pair=ks_pairwise_cor(Data=Ara,log=F)

# Reformate the output for upfollowing analysis
# Data is a list provided as an output of the function ks_pairwise_cor()
Ara_pair_max=ks_find_max_cor(Data=Ara_pair)

#move back to the main directory to write the output files
setwd("../")

##########################################################
# Create output files

# ks_make_table creates a data.frame with statistics per metabolite for the categories pairs, triplets and quadruples
# Please refer to the user's manual or manuscript for more details. 
# pair is a data.frame with the maximal correlated pairs generated by ks_find_max_cor
# Corr is a list with the maximal correlated triplets and quadruples generated by ks_find_max_cor_tr and ks_find_max_cor_qu
# Names are the names of the metabolites 
# tr is the correlation threshold, the pairs, triplets, and quadruples should satisfy the threshold to be considered in the output

write.table(ks_make_table(pair=Ecoli_pair_log_max,Corr=Ecoli_MaxCor,Names=rownames(Ecoli),tr=0.8),file = "Ecoli_complete_table_08.tab",row.names=F,col.names=T,quote = F,sep="\t")
write.table(ks_make_table(pair=Ecoli_pair_log_max,Corr=Ecoli_MaxCor,Names=rownames(Ecoli),tr=0.85),file = "Ecoli_complete_table_085.tab",row.names=F,col.names=T,quote = F,sep="\t")
write.table(ks_make_table(pair=Ecoli_pair_log_max,Corr=Ecoli_MaxCor,Names=rownames(Ecoli),tr=0.9),file = "Ecoli_complete_table_09.tab",row.names=F,col.names=T,quote = F,sep="\t")

write.table(ks_make_table(pair=Ara_pair_log_max,Corr=Ara_MaxCor,Names=rownames(Ara),tr=0.8),file = "Ara_complete_table_08.tab",row.names=F,col.names=T,quote = F,sep="\t")
write.table(ks_make_table(pair=Ara_pair_log_max,Corr=Ara_MaxCor,Name=rownames(Ara),tr=0.85),file = "Ara_complete_table_085.tab",row.names=F,col.names=T,quote = F,sep="\t")
write.table(ks_make_table(pair=Ara_pair_log_max,Corr=Ara_MaxCor,Names=rownames(Ara),tr=0.9),file = "Ara_complete_table_09.tab",row.names=F,col.names=T,quote = F,sep="\t")


############################################################    Comparative analysis
# Proceed from here if you want to compare the two species:

# ks_make_bipartite_graph creates a bipartite graph using the igraph package
# pairs is a data.frame with the maximal correlation of pairs generated by ks_find_max_cor
# triplets is a data.frame maximal correlated triplets generated by ks_find_max_cor_tr
# quadruples is a data.frame maximal correlated quadruples generated by ks_find_max_cor_qu
# tr is the correlation threshold, the pairs, triplets, and quadruples should satisfy the threshold to be considered in the output

Ecoli_bipartite_graph_08 = ks_make_bipartite_graph(pairs=Ecoli_pair_log_max,triplets=Ecoli_MaxCor$triplets,quadruples=Ecoli_MaxCor$quadruples,tr=0.8)
Ecoli_bipartite_graph_085 = ks_make_bipartite_graph(pairs=Ecoli_pair_log_max,triplets=Ecoli_MaxCor$triplets,quadruples=Ecoli_MaxCor$quadruples,tr=0.85)
Ecoli_bipartite_graph_09 = ks_make_bipartite_graph(pairs=Ecoli_pair_log_max,triplets=Ecoli_MaxCor$triplets,quadruples=Ecoli_MaxCor$quadruples,tr=0.9)

Ara_bipartite_graph_08 = ks_make_bipartite_graph(pairs=Ara_pair_log_max,triplets=Ara_MaxCor$triplets,quadruples=Ara_MaxCor$quadruples,tr=0.8)
Ara_bipartite_graph_085 = ks_make_bipartite_graph(pairs=Ara_pair_log_max,triplets=Ara_MaxCor$triplets,quadruples=Ara_MaxCor$quadruples,tr=0.85)
Ara_bipartite_graph_09 = ks_make_bipartite_graph(pairs=Ara_pair_log_max,triplets=Ara_MaxCor$triplets,quadruples=Ara_MaxCor$quadruples,tr=0.9)

# The function takes the bipartite graphs from the function ks_make_bipartite_graph from two different data set and creates an appropriate output,
# listing the metabolites and their corresponding coupling degree
# graph1: a bipartite graph of the first data set, generated by ks_make_bipartite_graph
# graph2: a bipartite graph of the second data set, generated by ks_make_bipartite_graph
# name1: the names of the metabolites of the first data set
# name2: the names of the metabolites of the second data set
# column_name1: string which is used to change the column name to match the data set
# column_name2: string which is used to change the column name to match the data set

data_08 = ks_graph_to_dataframe(Ecoli_bipartite_graph_08,Ara_bipartite_graph_08,rownames(Ecoli),rownames(Ara),"Ecoli metabolite coupling degree","Ara metabolite coupling degree")

data_085 = ks_graph_to_dataframe(Ecoli_bipartite_graph_085,Ara_bipartite_graph_085,rownames(Ecoli),rownames(Ara),"Ecoli metabolite coupling degree","Ara metabolite coupling degree")

data_09 = ks_graph_to_dataframe(Ecoli_bipartite_graph_09,Ara_bipartite_graph_09,rownames(Ecoli),rownames(Ara),"Ecoli metabolite coupling degree","Ara metabolite coupling degree")

out_list = list(data_08,data_085,data_09)

names(out_list) <- c("Comparison_08","Comparison_085","Comparison_09")
write_list(out_list,wb_name = "Metabolite_coupling_degree.xlsx")

###############################

# ks_shared_metabolites creates a list presenting the overlapping pairs, triplets and quadruples at the choosen threshold
# data_pair is a data.frame with the maximal correlation of pairs generated by ks_find_max_cor of the first data set
# data_all is a list with the maximal correlated triplets and quadruples generated by ks_find_max_cor_tr and ks_find_max_cor_qu of the firs data set
# data2_pair is a data.frame with the maximal correlation of pairs generated by ks_find_max_cor of the second data set
# data2_all is a list with the maximal correlated triplets and quadruples generated by ks_find_max_cor_tr and ks_find_max_cor_qu of the second data set
# tr is the correlation threshold, the pairs, triplets and quadruples should have at least to be considered for the output
# name1 is the string used to name the output of the first data set
# name2 is the string used to name the output of the second data set

Shared_metabolites_08 = ks_shared_metabolites(data_pair=Ecoli_pair_log_max,data_all=Ecoli_MaxCor,data2_pair=Ara_pair_log_max,data2_all=Ara_MaxCor,
                                              tr=0.8,name1="Ecoli",name2="Ara")
Shared_metabolites_085 = ks_shared_metabolites(data_pair=Ecoli_pair_log_max,data_all=Ecoli_MaxCor,data2_pair=Ara_pair_log_max,data2_all=Ara_MaxCor,
                                               tr=0.85,name1="Ecoli",name2="Ara")
Shared_metabolites_09 = ks_shared_metabolites(data_pair=Ecoli_pair_log_max,data_all=Ecoli_MaxCor,data2_pair=Ara_pair_log_max,data2_all=Ara_MaxCor,
                                              tr=0.9,name1="Ecoli",name2="Ara")

write_list(Shared_metabolites_08,wb_name ="Shared_metabolites_08.xls")
write_list(Shared_metabolites_085,wb_name ="Shared_metabolites_085.xls")
write_list(Shared_metabolites_09,wb_name ="Shared_metabolites_09.xls")


################################
# Save the results of the analysis
save.image("Run.RData")

# Clean up the directory if desired
unlink("Ara_run/", recursive = T)
unlink("Ecoli_run/", recursive = T)

