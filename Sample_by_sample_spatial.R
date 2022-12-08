"Script5: Sample by sample: spatial data comparison"
 

library('spatialLIBD')
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggridges)


#Load downloaded R object with all spatial data. Gene expression and layer annotations.
#http://research.libd.org/spatialLIBD/
sce <- fetch_data(type = 'sce')
#Get the log count matrix with all cells
log_count_matrix <- data.matrix(sce@assays@data@listData[["logcounts"]])
#Remove the no annotated cells
no_annotated <- metadata %>% filter(layer_guess!="NA" ) %>% select(barcode) 
rownames(no_annotated) <- NULL
no_annotated <- no_annotated$barcode
reduced_matrix <- log_count_matrix[,no_annotated]

#Identify the LR from previous analysis
all_LR_interactions <- read.table("/Users/mastorkia/Dropbox (EinsteinMed)/scRNASeq_Projects/AUTISM/CellChat/CORTEX/PRECORTEX_SAMPLE_BY_SAMPLE/ALL_ALL_interaction_results.csv", header = T)

gene_metadata <- data.frame(sce@rowRanges@elementMetadata@listData)
list_all_LR_interactions <-all_LR_interactions$pathways
#Using their data calculate the expression correlation of sygnificant L-R through all the cells (all matrix) or cell type wise. 
res_each_interaction <- data.frame()
for (i in list_all_LR_interactions) {
  test_lilist_all_LR_interactions <- data.frame(i)
  separated_lr <- separate(test_lilist_all_LR_interactions, 1, c("ligand","receptor"), sep="_")
  if(isFALSE(separated_lr$ligand==separated_lr$receptor)) {
    #Now using the ligand and receptor I mapp them in the gene_megadata to get the gene_id 
    desired_ligand <- gene_metadata %>% filter(gene_name==separated_lr$ligand) %>% select(gene_id)
    if(isTRUE(nrow(desired_ligand)>=1)) {
      ligand_geneid <- desired_ligand$gene_id
      count_ligand <- reduced_matrix[ligand_geneid,]
      desired_receptor <- gene_metadata %>% filter(gene_name==separated_lr$receptor) %>% select(gene_id)
      if(isTRUE(nrow(desired_receptor)>=1)){
        receptor_geneid <- desired_receptor$gene_id
        count_receptor <- reduced_matrix[receptor_geneid,]
        res <- cor.test(count_ligand, count_receptor, method = "pearson")
        group_df <- data.frame("Significant")
        combinated_df <- cbind(group_df,res$estimate)
        rownames(combinated_df) <- NULL
        combinated_df <- combinated_df %>% dplyr::rename(Group=1)%>% dplyr::rename(Correlation=2)
        res_each_interaction <- rbind(res_each_interaction,combinated_df)
      }
    }
  }
}


#Plot the density plot
res_each_interaction[is.na(res_each_interaction)] <- 0
density_plot <- ggplot(res_each_interaction, aes(x=Correlation)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") 
density_plot + theme(text = element_text(size = 14))  

#Plot the histogram but with cell by cell
layer_res_each_interaction <- data.frame()
cell_types <- unique(metadata$layer_guess)
cell_types <- cell_types[1:7]
for (i in cell_types) {
  layer_cells <- metadata %>% filter(layer_guess==i) %>% select(barcode) 
  rownames(layer_cells) <- NULL
  layer_cells <- layer_cells$barcode
  layer_cells_matrix <- reduced_matrix[,layer_cells]
  sig_each_interaction <- data.frame()
  for (j in list_all_LR_interactions) {
    test_lilist_all_LR_interactions <- data.frame(j)
    separated_lr <- separate(test_lilist_all_LR_interactions, 1, c("ligand","receptor"), sep="_")
    if(isFALSE(separated_lr$ligand==separated_lr$receptor)) {
      #Now using the ligand and receptor mapp them in the gene_megadata to get the gene_id 
      desired_ligand <- gene_metadata %>% filter(gene_name==separated_lr$ligand) %>% select(gene_id)
      if(isTRUE(nrow(desired_ligand)>=1)) {
        ligand_geneid <- desired_ligand$gene_id
        count_ligand <- layer_cells_matrix[ligand_geneid,]
        desired_receptor <- gene_metadata %>% filter(gene_name==separated_lr$receptor) %>% select(gene_id)
        if(isTRUE(nrow(desired_receptor)>=1)){
          receptor_geneid <- desired_receptor$gene_id
          count_receptor <- layer_cells_matrix[receptor_geneid,]
          res <- cor.test(count_ligand, count_receptor, method = "pearson")
          layer_df <- data.frame(paste0(i))
          group_df <- data.frame("Significant")
          combinated_df <- cbind(group_df,layer_df)
          first_df <- cbind(combinated_df,res$estimate)
          rownames(first_df) <- NULL
          first_df <- first_df %>% dplyr::rename(Group=1) %>% dplyr::rename(Layer=2) %>% dplyr::rename(Correlation=3)
          first_df <- subset(first_df, Correlation!="NA" )
          sig_each_interaction <- rbind(sig_each_interaction,first_df)
        }
      }
    }
  }
  layer_res_each_interaction <- rbind(layer_res_each_interaction,sig_each_interaction)
}

layer_res_each_interaction[is.na(layer_res_each_interaction)] <- 0
density_plot_by_groups <- ggplot(layer_res_each_interaction, 
                                 aes(x=Correlation, colour = Layer)) +  
  # Histogram with density instead of count on y-axis
  geom_density() 
density_plot_by_groups

#Now see if for randonmly selected genes the expression correlation keeps similar to significant L-R or not. To randomize the selection of sample use sample package
random_gene_list_1_1 <- data.frame(sample(gene_id, size=6000,replace = TRUE))
random_gene_list_1_1 <- random_gene_list_1_1 %>% dplyr::rename(gene_id=1)
random_gene_list_1_2 <- data.frame(sample(gene_id, size=6000,replace = TRUE))
random_gene_list_1_2 <- random_gene_list_1_2 %>% dplyr::rename(gene_id=1)
random_gene_list_1_3 <- data.frame(sample(gene_id, size=6000,replace = TRUE))
random_gene_list_1_3 <- random_gene_list_1_3 %>% dplyr::rename(gene_id=1)
random_gene_list <- rbind(random_gene_list_1_1,random_gene_list_1_2)
random_gene_list <- rbind(random_gene_list,random_gene_list_1_3)
random_gene_list_2 <- random_gene_list$gene_id
random_ligand_receptor <- data.frame() 
while(nrow(random_ligand_receptor) < 251) {
  gene1 <- sample(random_gene_list_2, size=1,replace = TRUE)
  gene2 <- sample(random_gene_list_2, size=1,replace = TRUE)
  random_l_r <- paste(gene1, gene2, sep ="_")
  random_ligand_receptor <- rbind(random_ligand_receptor,random_l_r )
  random_ligand_receptor <- random_ligand_receptor %>% dplyr::rename("Ligand_Receptor"=1)
}

#Loop through alll matrix and layers.
random_res_each_interaction <- data.frame()
random_ligand_receptor2 <- random_ligand_receptor$Ligand_Receptor
for (i in random_ligand_receptor2) {
  test_lilist_all_LR_interactions <- data.frame(i)
  separated_lr <- separate(test_lilist_all_LR_interactions, 1, c("ligand","receptor"), sep="_")
  if(isFALSE(separated_lr$ligand==separated_lr$receptor)) {
    #Using the ligand and receptor mapp them in the gene_megadata to get the gene_id 
    desired_ligand <- data.frame(separated_lr$ligand)
    if(isTRUE(nrow(desired_ligand)>=1)) {
      ligand_geneid <- desired_ligand$separated_lr.ligand
      count_ligand <- reduced_matrix[ligand_geneid,]
      desired_receptor <- data.frame(separated_lr$receptor)
      if(isTRUE(nrow(desired_receptor)>=1)){
        receptor_geneid <- desired_receptor$separated_lr.receptor
        count_receptor <- reduced_matrix[receptor_geneid,]
        res <- cor.test(count_ligand, count_receptor, method = "pearson")
        group_df <- data.frame("Random")
        combinated_df <- cbind(group_df,res$estimate)
        rownames(combinated_df) <- NULL
        combinated_df <- combinated_df %>% dplyr::rename(Group=1)%>% dplyr::rename(Correlation=2)
        random_res_each_interaction <- rbind(random_res_each_interaction,combinated_df)
      }
    }
  }
}
sum(is.na(random_res_each_interaction$Correlation))


#Plot the density plot
random_res_each_interaction[is.na(random_res_each_interaction)] <- 0
random_density_plot <- ggplot(random_res_each_interaction, 
                              aes(x=Correlation)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", 
                 fill="white") + geom_density(alpha=.2, fill="#FF6666") 
random_density_plot


#Now Plot the histogram but with cell by cell
random_layer_res_each_interaction <- data.frame()
for (i in cell_types) {
  layer_cells <- metadata %>% filter(layer_guess==i) %>% select(barcode) 
  rownames(layer_cells) <- NULL
  layer_cells <- layer_cells$barcode
  layer_cells_matrix <- reduced_matrix[,layer_cells]
  random_sig_each_interaction <- data.frame()
  for (j in random_ligand_receptor2) {
    test_lilist_all_LR_interactions <- data.frame(j)
    separated_lr <- separate(test_lilist_all_LR_interactions, 1, c("ligand","receptor"), sep="_")
    if(isFALSE(separated_lr$ligand==separated_lr$receptor)) {
      desired_ligand <- data.frame(separated_lr$ligand)
      if(isTRUE(nrow(desired_ligand)>=1)) {
        ligand_geneid <- desired_ligand$separated_lr.ligand
        count_ligand <- layer_cells_matrix[ligand_geneid,]
        desired_receptor <- data.frame(separated_lr$receptor)
        if(isTRUE(nrow(desired_receptor)>=1)){
          receptor_geneid <- desired_receptor$separated_lr.receptor
          count_receptor <- layer_cells_matrix[receptor_geneid,]
          res <- cor.test(count_ligand, count_receptor, method = "pearson")
          layer_df <- data.frame(paste0(i))
          group_df <- data.frame("Random")
          combinated_df <- cbind(group_df,layer_df)
          first_df <- cbind(combinated_df,res$estimate)
          rownames(first_df) <- NULL
          first_df <- first_df %>% dplyr::rename(Group=1) %>% dplyr::rename(Layer=2) %>% dplyr::rename(Correlation=3)
          first_df <- subset(first_df, Correlation!="NA" )
          random_sig_each_interaction <- rbind(random_sig_each_interaction,first_df)
        }
      }
    }
  }
  random_layer_res_each_interaction <- rbind(random_layer_res_each_interaction,random_sig_each_interaction)
}
random_layer_res_each_interaction[is.na(random_layer_res_each_interaction)] <- 0
random_density_plot_by_groups <- ggplot(random_layer_res_each_interaction, aes(x=Correlation, colour = Layer)) +  
  # Histogram with density instead of count on y-axis
  geom_density() 
random_density_plot_by_groups

#Plot random and significant densities together for each layer
density_results <- rbind(res_each_interaction,random_res_each_interaction)
all_histogram <- ggplot(density_results, aes(x=Correlation, colour = Group, fill=Group)) +
  geom_density(alpha=0.4)
all_histogram + theme(text = element_text(size = 14))  


#Plot random and significant densities together for each layer
density_results_all <- rbind(layer_res_each_interaction,random_layer_res_each_interaction)
all_histogram_layers <- ggplot(density_results_all) +
  geom_density_ridges(aes(x = Correlation, y = Layer,
                          group = interaction(Layer, Group),
                          fill = Group), scale = 1, alpha=0.4)
all_histogram_layers + theme(text = element_text(size = 14))  
