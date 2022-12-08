"Script3: Sample by sample approach: Ligand-Receptor interactions"
#Sample by sample approach. Finding dysregulated L-R interactions#
  
library(stringr)
library(dplyr)
library(plyr)
library(data.table)
library(fuzzySim)
library(graphics)
library(rstatix)
library(pheatmap)

#Loop to identify significant L-R interactions for desired pathways#
#Call files obtained in the global analysis.
my_files <- list.files("/Users/mastorkia/Dropbox (EinsteinMed)/scRNASeq_Projects/AUTISM/CellChat/CORTEX/PRECORTEX_SAMPLE_BY_SAMPLE/ALL_COMMUNICATIONS_LR", pattern="net")
my_files <- str_extract(my_files, '.*(?=\\.csv)')
#All the files in same format to read it correctly
for (i in my_files) {
  filename<-paste0("L_R_",i)
  wd<-paste0("ALL_COMMUNICATIONS_LR/",i,".csv")
  assign(filename,read.csv(wd, header=T, sep="\t"))
}


#Change the names of prob by each samples name.
lst <-  mget(ls(pattern="L_R_"))
lst1 <- Map(function(x,y) {names(x)[5] <- y; x},  lst, names(lst))
list2env(lst1, envir=.GlobalEnv)


total_big_merged_output <- data.frame()
new_name_files <- mget(ls(pattern="L_R_"))
for (i in new_name_files) {
  try2 <- i %>% select(1,2,7,5,9)
  total_big_merged_output <- join(total_big_merged_output,try2,type="full")
}


#Run Wilcoxon test and keep only interactions p < 0.05.
big_merged_output <- total_big_merged_output %>% select(!5)
source <- unique(big_merged_output$source)
target <- unique(big_merged_output$target)
interaction_results <- data.frame()
fdr_interaction_results <- data.frame()
for (i in source){
  for(j in target){
    try1 <- big_merged_output %>% filter(source==i & target==j) %>% select(3:26)
    rownames(try1) <- try1$interaction_name
    try1$interaction_name <- NULL
    transpose_try1 <- t(try1)
    transpose_try1 <- data.frame(transpose_try1)
    Condition <- c("ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control")
    try2<-cbind(Condition,transpose_try1)
    rownames(try2) <- NULL
    try2<-try2[, which(colMeans(!is.na(try2)) > 0.32)]
    if (isTRUE(ncol(try2) >= 2)) {interaction_output_ttest <- list()
    interaction_ttest<-colnames(try2)
    interaction_ttest <-interaction_ttest[-1]
    #This one only works for wilcox_test p values
    for(x in interaction_ttest){
      fmla <- formula(paste0(x, " ~ Condition"))
      res<-wilcox.test(fmla, data = try2, exact = FALSE)
      interaction_output_ttest[[x]] <- res$p.value
    }
    interaction_output_ttest <- unlist(interaction_output_ttest)
    interaction_output_ttest<-data.frame(interaction_output_ttest)
    interaction_output_ttest$pathways<-rownames(interaction_output_ttest)
    rownames(interaction_output_ttest)<- NULL
    interaction_output_ttest <- interaction_output_ttest %>% 
      dplyr::rename(!!paste0(i,"_",j):=interaction_output_ttest)
    #filter(interaction_output_ttest < 0.05) %>%
    interaction_results <- join(interaction_results,interaction_output_ttest,type="full")
    }
  }
}


#Prepare it to plot it in a heatmap
all_na <- function(x) any(!is.na(x))
interaction_results2 <- interaction_results %>% select_if(all_na)
rownames(interaction_results2) <- interaction_results2$pathways
interaction_results2$pathways <- NULL
interaction_results2_matrix <- data.matrix(interaction_results2)
ordered.try1 <- interaction_results2_matrix[order(rownames(interaction_results2_matrix)),]   #this orders the rownames
ordered.try1 <- ordered.try1[,order(colnames(ordered.try1))]   #this orders the by colnames
colnames(ordered.try1) <- gsub("Oligodendrocytes", "Oligo", colnames(ordered.try1))
sig_cell_cell_heatmap = pheatmap(ordered.try1, color = hcl.colors(50, "BluYl"),cluster_rows = F, cluster_cols = F, fontsize_row = 5, fontsize_col = 5, na_col="white" )
sig_cell_cell_heatmap

#Heatmap showing which L-R have a negative or positive fold change.
heatmap_ordered.try1 <- ordered.try1
for (i in rownames(heatmap_ordered.try1)){
  see <-heatmap_ordered.try1[i,] 
  see2 <- names(na.omit(see))
  for(j in see2){
    cell_names <- str_split(j,"_", simplify = T)
    interaction_strengh <- big_merged_output %>% 
      filter(interaction_name==i) %>% 
      filter(source==cell_names[1,1]) %>%
      filter(target==cell_names[1,2])
    fold_change <- rowMeans(interaction_strengh %>% select(contains("asd")), na.rm = T) - rowMeans(interaction_strengh %>% select(contains("Control")), na.rm = T)
    ifelse(fold_change>0, heatmap_ordered.try1[i,j] <- -(log(heatmap_ordered.try1[i,j])), heatmap_ordered.try1[i,j] <- log(heatmap_ordered.try1[i,j]))
  }
}

#Plot again the heatmap
fold_change_sig_cell_cell_heatmap = pheatmap(heatmap_ordered.try1, color = colorRampPalette(c("blue", "white", "red"))(100),cluster_rows = F, cluster_cols = F, fontsize_row = 5, fontsize_col = 5, na_col="white" )
fold_change_sig_cell_cell_heatmap

#Number of rows with up or down interactions and how many cell pairs are involved in each L-R interaction
up_down_results <- data.frame()
for (i in rownames(heatmap_ordered.try1)){
  ole <-heatmap_ordered.try1[i,] 
  ole2 <- na.omit(ole)
  hi <- data.frame(length(ole2))
  hi2 <- sum(ole2)
  hi3 <- cbind(hi, hi2)
  rownames(hi3) <- i
  up_down_results <- rbind(up_down_results, hi3)
}
up_down_results <- up_down_results %>% dplyr::rename(n_cell_cell_pairs=length.ole2.) %>% dplyr::rename(overall_strength=hi2)
table(up_down_results$overall_strength > 0 )["TRUE"]
table(up_down_results$overall_strength > 0 )["FALSE"]


#Heatmap having the upregulated interactions in ASD all together on the top and separated by pathways
pathways_lr <- total_big_merged_output %>% select(pathway_name, interaction_name) %>% unique()
heatmap_ordered.try2 <- as.data.frame(heatmap_ordered.try1)
heatmap_ordered.try2$interaction_name <- rownames(heatmap_ordered.try2)
rownames(heatmap_ordered.try2) <- NULL
sig_path_lr <- inner_join(heatmap_ordered.try2,pathways_lr, type="inner" )

test <- sig_path_lr %>% replace(is.na(.), 0) %>% mutate(row_wise_sum = rowSums(.[1:155])) %>% group_by(pathway_name) %>% summarise_at(vars(row_wise_sum),funs(sum(.,na.rm=TRUE))) %>% arrange(desc(row_wise_sum))
test
#Now order the data frame based on the results of test
si <- sig_path_lr[order(match(sig_path_lr$pathway_name, ordering)), ]

#Again create the data matrix for heatmap
sig_path_lr2 <- si
rownames(sig_path_lr2) <- sig_path_lr2$interaction_name
sig_path_lr2$pathway_name <- NULL
sig_path_lr2$interaction_name <- NULL
sig_path_lr3 <- data.matrix(sig_path_lr2)


#Now plot the heatmap
specific_ordered_fold_change_sig_cell_cell_heatmap = pheatmap(sig_path_lr3, color = colorRampPalette(c("blue", "white", "red"))(100),cluster_rows = F, cluster_cols = F, fontsize_row = 5, fontsize_col = 5, na_col="white" )
specific_ordered_fold_change_sig_cell_cell_heatmap


#With all the p values from interaction_results draw a p value histogram,
interaction_results_histogram <- interaction_results %>% select (!2)
interaction_results_histogram <- data.matrix(interaction_results_histogram)
h <- hist(interaction_results_histogram,
          main="all L-R p values",
          xlab="p values",
          freq=T, 
          ylab = "Freq", 
          ylim = c(0,900), 
          breaks=30
)
h


#FDR multiple testing by cell-cell pair families
interaction_results2 <- interaction_results
interaction_results2$pathways <- NULL
colnames(interaction_results2) <- gsub("/", "_", colnames(interaction_results2))
colnames(interaction_results2) <- gsub("-", "_", colnames(interaction_results2))
fdr_p_results <- data.frame(ncol = 1, nrow = 270)
pathways_list <- colnames(interaction_results2)
for (i in colnames(interaction_results2)){
  vector1 <- interaction_results2[[i]]
  results_p <-p.adjust(vector1, method = "fdr")
  results_p1 <- data.frame(results_p)
  if (isTRUE(nrow(results_p1) >= 1)) {
    results_p1 <- results_p1 %>%  dplyr::rename(!!paste0(i):=results_p)
    fdr_p_results <-cbind(fdr_p_results,results_p1)
  }
}

#fdr_p_results <- fdr_p_results %>% select(-(1:2))
#fdr_p_results_cutoff <- fdr_p_results %>% filter(fdr_p_results < 0.25)

#write.table(fdr_p_results, "fdr_results.csv", col.names = T, row.names = F, quote = F, sep = "\t")

#Extract all the networks and L_R through each pathway list.
pathway_list_complete <- read.table("PATHWAY_LIST.txt")
pathway_list_complete <- pathway_list_complete$V1
pathways_gene<-list()
for (i in pathway_list_complete){
  ligands <- data.frame()
  receptors <- data.frame()
  for(j in new_name_files){
    test<-j %>% dplyr::filter(pathway_name==i) 
    test <- unique(test$receptor)
    test <- data.frame(test)
    ligands <- join(ligands, test, type="full")
    test2<-j %>% dplyr::filter(pathway_name==i) 
    test2<- unique(test2$ligand) 
    test2 <- data.frame(test2)
    receptors <- join(receptors, test2, type="full")
  }
  ligands2 <- t(ligands)
  receptors2 <- t(receptors)
  all_l_r <- cbind(ligands2,receptors2)
  pathways_gene[[i]] <- all_l_r
}
table_pathways_gene <-unlist(pathways_gene)
write.table(table_pathways_gene, "table_pathways_gene.csv", col.names = F, row.names = T, quote = F, sep = "\t")


