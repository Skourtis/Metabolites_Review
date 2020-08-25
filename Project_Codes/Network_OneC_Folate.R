library(pacman)
p_load("plyr", "dplyr", "openxlsx","tidyverse", "Hmisc", "corrplot","igraph","RCy3","BioNet","fuzzyjoin", 
       "naniar", "plotly", "reshape", "KEGGgraph", "KEGG.db", "Rgraphviz", "RBGL", "KEGGprofile", "visNetwork", "ggnetwork",
       "rlist", "leaflet", "matrixStats", "readxl", "RColorBrewer", "pheatmap")


map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

Protein_Families <- as.data.frame(read_tsv("./Project_Datasets/uniprot-organism__Homo+sapiens+(Human)+[9606]_.tab")) 
Protein_Families$Gene_names2 <- str_match(Protein_Families$Gene_names, pattern = "^([:graph:]*) ")[,2] 
Protein_Families$Gene_names[!is.na(Protein_Families$Gene_names2)] <- na.omit(Protein_Families$Gene_names2)
Protein_Families$Protein_families2 <- str_match(Protein_Families$Protein_families, pattern = "^([[:print:]]*?),")[,2] 
Protein_Families$Protein_families[!is.na(Protein_Families$Protein_families2)] <- na.omit(Protein_Families$Protein_families2)
Protein_Families$Protein_families2 <- str_match(Protein_Families$Protein_families, pattern = "^([[:print:]]*?);")[,2] 
Protein_Families$Protein_families[!is.na(Protein_Families$Protein_families2)] <- na.omit(Protein_Families$Protein_families2)
KEGG_to_Gene_name <- openxlsx::read.xlsx("./Project_datasets/uniprot_hsa01100.xlsx")[,c(1,6)] %>%
    mutate(Gene.names = str_match(.$Gene.names, pattern = "^([:graph:]*)")[,2])  %>% na.omit()
KEGG_to_Gene_name <- KEGG_to_Gene_name[!duplicated(KEGG_to_Gene_name$KEGG_ID),]
###Pertubation score
Pertubation_score <- read.csv("./Project_Datasets/Pertubation_score.csv")
Pertubation_score <- data.frame(Gene_name =str_remove_all(Pertubation_score$X, "\\([:graph:]*\\)"),
                                     Color = rowQuantiles(as.matrix(Pertubation_score[,-1]), na.rm = T)[,2]) %>% 
    mutate(width = 1/rowSds(as.matrix(Pertubation_score[,-1]), na.rm = T))

### Retrieving and Merging KEGG Folate
tmp <- tempfile()
genes_reactions <- data.frame(Gene_id = NULL,
                              Reaction = NULL, 
                              pathway = NULL,
                              stringsAsFactors = F)
compound_reactions <- data.frame(Reaction = NULL,
                                 Substrate = NULL,
                                 Product = NULL,
                                 Direction = NULL,
                                 pathway = NULL,
                                 stringsAsFactors = F)


pathways <-  c("hsa00670", "hsa00240", "hsa00230","hsa00270")
for (j in pathways){
    retrieveKGML(j, organism="hsa", destfile=tmp, method="auto", quiet=TRUE)
    testing <- KEGGgraph::parseKGML2Graph(tmp)
    
    test_folate <- KEGGgraph::parseKGML(tmp)
    # mapkG2 <- KEGGpathway2Graph(test_folate, expandGenes=TRUE)
    # set.seed(124)
    # randomNodes <- sample(nodes(mapkG2), 25)
    # mapkGsub <- subGraph(randomNodes, mapkG2)
    # plot(mapkGsub)
    
    for (i in 1:length(testing@nodeData@defaults[["KEGGNode"]][["nodes"]])){
        df <- data.frame(Gene_id = testing@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@entryID,
                         Reaction = testing@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@reaction,
                         pathway = j,
                         stringsAsFactors = F)
        genes_reactions <- rbind(genes_reactions,df)
        
    }
    
    
    for (i in 1:length(test_folate@reactions)){
        df <- data.frame(Reaction = test_folate@reactions[[i]]@name,
                         Substrate = test_folate@reactions[[i]]@substrateName,
                         Product = test_folate@reactions[[i]]@productName,
                         Direction = test_folate@reactions[[i]]@type,
                         pathway = j,
                         stringsAsFactors = F)
        compound_reactions <- rbind(compound_reactions,df)
        
    }
    
    
}
Enzymes_Metabolites <- read_tsv("./Project_Datasets/all_unique_KEGG_metabolites_mapped_KEGG.tsv")
compound_reactions <- separate_rows(compound_reactions, Reaction, sep =" ") %>% 
    mutate(Substrate = str_remove(Substrate,"cpd:"),
           Product = str_remove(Product,"cpd:"))
genes_reactions <- separate_rows(genes_reactions, Reaction, sep =" ")
Pathways <- full_join(genes_reactions,compound_reactions, by = c("Reaction", "pathway")) %>% unique()
Pathways <- left_join(Pathways,KEGG_to_Gene_name, by= c("Gene_id" = "KEGG_ID")) %>% 
    left_join(distinct(Enzymes_Metabolites[,1:2]) , by = c( "Substrate"  ="KEGG" )) %>% 
    left_join(distinct(Enzymes_Metabolites[,1:2]), by = c( "Product"  ="KEGG" ))  %>%
    dplyr::rename(From = `Compound.x`,  To = `Compound.y`) %>% 
    distinct( Gene_id, pathway, Substrate, Product, Direction, Gene.names, .keep_all = T)

Pathways$Gene.names[is.na(Pathways$Gene.names)] <- Pathways$Gene_id[is.na(Pathways$Gene.names)]
Pathways$From[is.na(Pathways$From)] <- Pathways$Substrate[is.na(Pathways$From)]
Pathways$To[is.na(Pathways$To)] <- Pathways$Product[is.na(Pathways$To)]

M00034_nodes <-  c("C00073", "C00019", "C01137", "C00170" ,"C03089",  "C04188",  "C04582","C15650" , "C15651" ,"C15606" , "C01180")
Pathways <- anti_join(Pathways,inner_join(Pathways[Pathways$pathway == "hsa00270" & !(Pathways$Substrate %in% M00034_nodes),],
                   Pathways[Pathways$pathway == "hsa00270" & !(Pathways$Product %in% M00034_nodes),]))


#####
nodes <- c(Pathways$From,Pathways$To) %>% unique()
nodes <- data.frame(id = nodes, label = nodes)
#colnames(nodes) <- c('id', "label")
#nodes$label[is.na(nodes$label)] <- nodes$id[is.na(nodes$label)]
nodes$title <- paste0("<p><b>",nodes$id ,"</b><br></p>")



edges_enzymes <- Pathways[,c( "From","To", "Gene.names", "Direction")] %>% mutate(font.color = "green", font.size = 10)%>%
    dplyr::rename(label = `Gene.names`, from = From, to = To) %>% left_join(Pertubation_score, by = c("label" = "Gene_name"))



edges_enzymes$width[is.na(edges_enzymes$width)] <- min(edges_enzymes$width, na.rm =  T)
grey_edges <- which(is.na(edges_enzymes$Color))
edges_enzymes$Color[is.na(edges_enzymes$Color)] <- mean(edges_enzymes$Color, na.rm =  T)
edges_enzymes$label[is.na(edges_enzymes$label)] <- " "
edges_enzymes$title[edges_enzymes$label == " "] <- " "


#Import Groups

Metabol_clusters <- data.frame(Metabolite = NULL, group = NULL)
for(i in 1:3){ 
    Metabol_clusters <- rbind(Metabol_clusters,data.frame (Metabolite = read_tsv(paste0("./Project_Datasets/metabolite_module_",i,".txt"))$V1, group = i))
}

mypal <- colorRampPalette( c( "red","purple", "#0080ff" ) )(nrow(edges_enzymes))


#pal <- colorNumeric("RdBu", domain = (min(edges_enzymes$Color)-1):(max(edges_enzymes$Color)+1), reverse = T)
edges_enzymes$color <- map2color(edges_enzymes$Color,mypal)
edges_enzymes$color[grey_edges] <- "#a6a6a6"
edges_enzymes$arrows.middle <- if_else(edges_enzymes$Direction == "irreversible", T, F)
edges_enzymes$dashes <- if_else(edges_enzymes$label == " ", T, F)
edges_enzymes <- left_join(edges_enzymes, na.omit(distinct(Protein_Families[,c(5,8)])), by = c("label" = "Gene_names")) 
edges_enzymes <- edges_enzymes %>% left_join(na.omit(ddply(edges_enzymes, .(from, to, Protein_families), nrow)), by =c("from", "to", "Protein_families"))
edges_enzymes$Gene <- edges_enzymes$label
edges_enzymes$V1[is.na(edges_enzymes$V1)] <- 1
edges_enzymes$label[edges_enzymes$V1>2] <- edges_enzymes$Protein_families[edges_enzymes$V1>2]
edges_enzymes <- distinct(edges_enzymes)

edges_enzymes$color[edges_enzymes$label == edges_enzymes$Protein_families] <- "green"
edges_enzymes$hidden <- if_else(duplicated(edges_enzymes[,c("from", "to", "Protein_families")]) & edges_enzymes$color == "green",T,F)
edges_enzymes$arrows.middle[edges_enzymes$label == edges_enzymes$Protein_families] <- FALSE
edges_enzymes$title <- paste0(edges_enzymes$Gene, " :: Pertubation : ", edges_enzymes$Color)
edges_enzymes <- edges_enzymes %>% 
    dplyr::group_by(Protein_families, from, to) %>% 
    mutate(title = paste0(title, collapse = "<br>"))

Metabolites_mapped <- read.csv("./Project_Datasets/MetaboAnalyst_mapping.csv") %>%  na.omit() # left_join(Metabol_clusters, by = c("Query"= "Metabolite")) %>%
    #.[,c("KEGG","group")] %>%
nodes <- left_join(nodes, Metabolites_mapped, by = c("id" = "KEGG"))
nodes$color.background <- "#f8f8ff"
edges_enzymes$label <- as.character(edges_enzymes$label)
edges_enzymes <- edges_enzymes [!duplicated(edges_enzymes[1:3]),]
visNetwork(nodes,edges_enzymes, height = 1000, width = "100%") %>% 
    visIgraphLayout(smooth = T)%>% 
    visNodes(
        shape = "box",
        shadow = list(enabled = TRUE, size = 10)
    ) %>%
         visEdges(smooth = list(type = "dynamic"), color = list(highlight = "#C62F4B",opacity = 0.35, border = "black"), 
             font = list(align = "middle")) %>%
    visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
               selectedBy = "label") %>%
    visPhysics(stabilization = FALSE,solver = "forceAtlas2Based", 
               forceAtlas2Based = list(gravitationalConstant = -90))


#heatmap
#Metabolites vs. Tissues#
CellLineInputFile="./Project_Datasets/41591_2019_404_MOESM2_ESM.xlsx"
CellLineAnnoFile="./Project_Datasets/mod_tissue_origins.txt"
FolateMetaFile="./Project_Datasets/kegg_metabolites.txt"

cell_line_data=read_excel(CellLineInputFile,sheet=3)
mod_cell_line_data=data.frame(column_to_rownames(cell_line_data, var = "...1"),check.names=F)

cell_line_annotation=read.delim(file=CellLineAnnoFile,header=F)
tissue_annotation=str_replace(cell_line_annotation[,1],"HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","BLOOD")
pat_others=c("AUTONOMIC_GANGLIA|BILIARY_TRACT|BONE|KIDNEY|PLEURA|PROSTATE|SALIVARY_GLAND|SOFT_TISSUE|THYROID")
tissue_annotation=str_replace_all(tissue_annotation,pat_others,"OTHERS")

mod_tissue_annotation=data.frame(tissue_annotation)
row.names(mod_tissue_annotation) <- rownames(mod_cell_line_data)

#Colour palette for the tissue annotation#
n <- 15
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#ann_colors=sample(col_vector, n)
n <- 15
ann_colors=col_vector[50:(50+n-1)]
names(ann_colors) <- unique(mod_tissue_annotation$tissue_annotation)
ann_colors=list(tissue_annotation=ann_colors)

#folate_metabolites=c("glycine", "serine", "adenosine", "glutamate", "NAD", "NADP", "methionine", "betaine", "cystathionine", "choline", "dimethylglycine", "homocysteine", "5-adenosylhomocysteine", "putrescine", "sarcosine", "alpha-ketoglutarate" )
folate_metabolites=read.delim(file=FolateMetaFile,header=F, stringsAsFactors = F)
folate_cell_line_data=subset(mod_cell_line_data,select= colnames(mod_cell_line_data) %in% folate_metabolites[,1])
folate_cell_line_data_test <- folate_cell_line_data %>% mutate(Cell_line = rownames(folate_cell_line_data)) %>%
    pivot_longer(-Cell_line, names_to= "Metabolite", values_to = "Abundance") %>% mutate(Organ = str_match(.$Cell_line, "_([:graph:]*$)")[,2])


#pheatmap(t(folate_cell_line_data),color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(50),scale="row",clustering_distance_rows="correlation",clustering_distance_cols = "correlation",show_rownames =T,show_colnames =F,annotation_col=mod_tissue_annotation,annotation_colors=ann_colors,main="Metabolites and Tissues in the Folate Pathway")

heatmap <- pheatmap(folate_cell_line_data,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(50),scale="column",clustering_distance_rows="correlation",clustering_distance_cols = "correlation",show_rownames =F,show_colnames =T,annotation_row=mod_tissue_annotation,annotation_colors=ann_colors,main="Metabolites and Tissues in the Folate Pathway")

save(nodes, #done
     edges_enzymes, #done
     folate_cell_line_data_test, #done
     Metabolites_mapped,#done
     Enzymes_Metabolites, #done
     heatmap, #done
     file = "./Project_Output/Trial_network.RData")
