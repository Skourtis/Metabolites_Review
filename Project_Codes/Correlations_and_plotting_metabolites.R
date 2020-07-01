library(pacman)
install.packages("remotes")
remotes::install_github("larmarange/JLutils")
library(larm)

p_load("openxlsx","tidyverse", "Hmisc", "corrplot","igraph","RCy3","BioNet")
library(BiocManager)
#BiocManager::install("RCy3")
Reactome_pathways <- read_tsv("./Project_Datasets/ChEBI2Reactome_PE_Pathway.txt", col_names = F) %>% 
    .[.$X8 == "Homo sapiens",]

cytoscapePing ()
cytoscapeVersionInfo ()
Metabolites <- read.xlsx("./Project_Datasets/CCLE_metabolites_landscape_of_cancer.xlsx",
                         sheet = "1-clean data")
df<-scale(Metabolites[,-1])# normalize the data frame. This will also convert the df to a matrix.  
Metabolites <- as.data.frame(t(Metabolites), stringsAsFactors = F) %>% mutate(Names = rownames(.))

dist_mat <- dist(t(df), method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)
cut_avg  <- cutree(hclust_avg, 10)
clustered_metabolites <- data.frame(Metabolite = names(cut_avg), Cluster = cut_avg)
clustered_metabolites<- clustered_metabolites%>% mutate(MSEA = case_when(
    Cluster == 2 ~ 'Glutamate Metabolism',
    Cluster == 1 ~ 'Tryptophan Metabolism',
    Cluster == 3 ~ 'Pyrimidine Metabolism',
    Cluster == 5 ~ 'Methionine Metabolism',
    TRUE ~ "Other"
)
)
clustered_metabolites <- left_join(clustered_metabolites,Metabolites, by = c("Metabolite" = "Names"))
clustered_metabolites$`929` <- as.numeric(clustered_metabolites$`1`)
clustered_metabolites$`930` <- as.numeric(clustered_metabolites$`2`)
clustered_metabolites$`931` <- as.numeric(clustered_metabolites$`52`)
loadTableData(
    clustered_metabolites[,c("Metabolite","931")],
    data.key.column = "Metabolite",
    table = "node",
    table.key.column = "name",
)

write.csv(clustered_metabolites, "./Project_Output/Clustered_metabolites_CCLE_for_MetaboAnalyst.csv")
test <- t(Metabolites)
corr<-rcorr(df) # compute Pearson's (or spearman's corr) with rcorr from Hmisc package. I like rcorr as it allows to separately access the correlations, the # or observations and the p-value. ?rcorr is worth a read.
corr_r<-as.matrix(corr[[1]])# Access the correlation matrix. 
#corr_r[,1]# subset the correlation of "a" (=var1 ) with the rest if you want.
pval<-as.matrix(corr[[3]])# get the p-values

corr_r[row(corr_r) == col(corr_r) ] <- 0

# set all correlations that are less than 0.9 to zero
corr_r[which(corr_r<0.5)] <- 0

#get rid of rows and columns that have no correlations with the above thresholds
corr_r <- corr_r[which(rowSums(corr_r) != 0),which(colSums(corr_r) !=0)]

#write out the correlation file
correlation_filename <- "./Project_Output/Metabolites_Matrix.txt"
write.table( corr_r,  file = correlation_filename, col.names  = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)

amat_url <- "aMatReader/v1/import"
amat_params = list(files = list("C:/Users/skourtis/OneDrive - CRG - Centre de Regulacio Genomica/Bioinformatics Projects/Metabolites_Review/Project_Output/Metabolites_Matrix.txt"),
                   delimiter = "TAB",
                   undirected = FALSE,
                   ignoreZeros = TRUE,
                   interactionName = "correlated with",
                   rowNames = FALSE
)

response <- cyrestPOST(operation = amat_url, body = amat_params, base.url = "http://localhost:1234")

current_network_id <- response$data["suid"]
layoutNetwork('cose',
              network = as.numeric(current_network_id))
renameNetwork(title = "Metabolites_Matrix",
              network = as.numeric(current_network_id))

nodes_in_network <- rownames(corr_r)

#make sure it is set to the right network
setCurrentNetwork(network = getNetworkName(suid=as.numeric(current_network_id)))

#cluster the network
clustermaker_url <- paste("cluster mcl network=SUID:",current_network_id, sep="")
commandsGET(clustermaker_url)

#get the clustering results
default_node_table <- getTableColumns(table= "node",network = as.numeric(current_network_id))

head(default_node_table)

#add an additional column to the gene scores table to indicate in which samples
# the gene is significant


# installation_responses <- c()
# 
# #list of app to install
# cyto_app_toinstall <- c("clustermaker2", "enrichmentmap", "autoannotate", "wordcloud", "stringapp", "aMatReader")
# 
# cytoscape_version <- unlist(strsplit( cytoscapeVersionInfo()['cytoscapeVersion'],split = "\\."))
# if(length(cytoscape_version) == 3 && as.numeric(cytoscape_version[1]>=3) 
#    && as.numeric(cytoscape_version[2]>=7)){
#     for(i in 1:length(cyto_app_toinstall)){
#         #check to see if the app is installed.  Only install it if it hasn't been installed
#         if(!grep(commandsGET(paste("apps status app=\"", cyto_app_toinstall[i],"\"", sep="")), 
#                  pattern = "status: Installed")){
#             installation_response <-commandsGET(paste("apps install app=\"", 
#                                                       cyto_app_toinstall[i],"\"", sep=""))
#             installation_responses <- c(installation_responses,installation_response)
#         } else{
#             installation_responses <- c(installation_responses,"already installed")
#         }
#     }
#     installation_summary <- data.frame(name = cyto_app_toinstall, 
#                                        status = installation_responses)
#     
#     knitr::kable(list(installation_summary),
#                  booktabs = TRUE, caption = 'A Summary of automated app installation'
#     )
# }


####MetaboAnalystR
##Rtools needed to be updated, and installed manually, and added to path with : 
#writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron") after R was restarted
packageurl <- "http://cran.r-project.org/src/contrib/Archive/mnormt/mnormt_1.5-7.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
pacman::p_load("crmn","impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea")
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)
library(MetaboAnalystR)

### Enrichment Analysis
tmp.vec <- c("Acetoacetic acid", "Beta-Alanine", "Creatine", "Dimethylglycine", "Fumaric acid", "Glycine", "Homocysteine", "L-Cysteine", "L-Isolucine", "L-Phenylalanine", "L-Serine", "L-Threonine", "L-Tyrosine", "L-Valine", "Phenylpyruvic acid", "Propionic acid", "Pyruvic acid", "Sarcosine")
# Create mSetObj
mSet<-InitDataObjects("conc", "msetora", FALSE)
#Set up mSetObj with the list of compounds
mSet<-Setup.MapData(mSet, tmp.vec);
# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
# mSet<-CrossReferencing(mSet, "name");
# 
# 
# # Example compound name map
# mSet$name.map
# $query.vec
# [1] "Acetoacetic acid" "Beta-Alanine" "Creatine" "Dimethylglycine" "Fumaric acid"
# [6] "Glycine" "Homocysteine" "L-Cysteine" "L-Isolucine" "L-Phenylalanine"
# [11] "L-Serine" "L-Threonine" "L-Tyrosine" "L-Valine" "Phenylpyruvic acid"
# [16] "Propionic acid" "Pyruvic acid" "Sarcosine"
# $hit.inx
# [1] 42 40 46 62 88 78 588 446 NA 104 120 109 103 702 131 159 164 185
# $hit.values
# [1] "Acetoacetic acid" "Beta-Alanine" "Creatine" "Dimethylglycine" "Fumaric acid"
# [6] "Glycine" "Homocysteine" "L-Cysteine" NA "L-Phenylalanine"
# [11] "L-Serine" "L-Threonine" "L-Tyrosine" "L-Valine" "Phenylpyruvic acid"
# [16] "Propionic acid" "Pyruvic acid" "Sarcosine"
# $match.state
# [1] 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1
# Continute with the enrichment analysis. . .
# # Create the mapping results table
# mSet<-CreateMappingResultTable(mSet)
# # Input the name of the compound without any matches
# mSet<-PerformDetailMatch(mSet, "L-Isolucine");
# # Create list of candidates to replace the compound
# mSet <- GetCandidateList(mSet);
# # Identify the name of the compound to replace
# mSet<-SetCandidate(mSet, "L-Isolucine", "L-Isoleucine");
# # Set the metabolite filter
# mSet<-SetMetabolomeFilter(mSet, F);
# # Select metabolite set library
# mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
# # Calculate hypergeometric score, results table generated in your working directory
# mSet<-CalculateHyperScore(mSet)
# # Plot the ORA, bar-graph
# mSet<-PlotORA(mSet, "ora_0_", "bar", "png", 72, width=NA)

exportNetwork(
    filename = "./Project_Output/Pathways_MSEA.tif",
    type = "tif"
)
