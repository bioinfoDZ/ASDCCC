"Script2: Sample by sample approach: network analysis"

#Sample by sample approach to find dysregulated networks#
  
#Loop to find dysregulated networks through all cell-cell interactions between samples. Strength of the interactions was calculated using CellChat for each sample. 
  
#Find the dysregulated networks between ASD and Control for each pair of all cell types.

library(stringr)
library(dplyr)
library(plyr)
library(data.table)
library(ggbiplot)
library("msos")
library(circlize)
library(reshape)
library(graph4lg)
library(Matrix)


#Load excel files with network probabilities for each of the samples that were obtained by CellChat.
my_files <- list.files(pattern="netP")
my_files <- str_extract(my_files, '.*(?=\\.csv)')
#All the files in same format to read it correctly
for (i in my_files) {
  filename<-paste0(i)
  wd<-paste0(i,".csv")
  assign(filename,read.csv(wd, header=T, sep="\t"))
}
lst_files <-  mget(ls(pattern="netP"))
lst1_files <- Map(function(x,y) {names(x)[4] <- y; x},  lst_files, names(lst_files))
list2env(lst1_files, envir=.GlobalEnv)

#Create a unique big file for all networks and all samples.
reduced_files <- data.frame()
lst_files <-  mget(ls(pattern="netP"))
for (i in lst_files) {
  try2 <- i %>% dplyr::select(1:4)
  reduced_files <- join(reduced_files,try2,type="full")
}
write.table(reduced_files, "probability_pathway_celltype.csv", quote = F, row.names = F, col.names = T, sep = "\t")


#Get the overall network wilcoxon results between ASD and Control samples. All cell communications together.
summed_pathways_big<-data.frame()
pathway_output<-data.frame()
for (i in lst_files){
  pathways<-unique(i$pathway_name)
  name<-colnames(i)[4]
  for(x in name){
    pathway_sum_output <- list()
    for(j in pathways){
      pathway_sum_output[[j]]<-sum(i[i$pathway_name==j,][,4])}
    pathway_sum_output <- unlist(pathway_sum_output)
    pathway_sum_output<-data.frame(pathway_sum_output)
    pathway_sum_output$pathways<-rownames(pathway_sum_output)
    rownames(pathway_sum_output)<- NULL
    pathway_sum_output <- pathway_sum_output %>% dplyr::rename(!!paste0(x):=pathway_sum_output)
    pathway_sum_output<- pathway_sum_output[,c(2,1)]
    summed_pathways_big<- join(summed_pathways_big,pathway_sum_output, type="full")
  }
}

#From here get the p values of Wilcoxon test for all networks ASD vs Control
rownames(summed_pathways_big) <- summed_pathways_big$pathways
summed_pathways_big$pathways <- NULL
transpose_summed_pathways_big <- t(summed_pathways_big)
transpose_summed_pathways_big <- data.frame(transpose_summed_pathways_big)
Condition <- c("ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control")
try2<-cbind(Condition,transpose_summed_pathways_big)
rownames(try2) <- NULL
try2<-try2[, which(colMeans(!is.na(try2)) > 0.3)]
if (isTRUE(ncol(try2) >= 2)) {
  cell_cell_output_ttest <- list()
  pathways_ttest<-colnames(try2)
  pathways_ttest <- pathways_ttest[-1]
  pathway_output_ttest <- data.frame()
  pathway_output_ttest <- list()
  for(x in pathways_ttest){
    fmla <- formula(paste0(x, " ~ Condition"))
    res<-wilcox.test(fmla, data = try2, exact = FALSE)
    pathway_output_ttest[[x]] <- res$p.value
  }
  pathway_output_ttest <- unlist(pathway_output_ttest)
  pathway_output_ttest<-data.frame(pathway_output_ttest)
  pathway_output_ttest$pathways<-rownames(pathway_output_ttest)
  rownames(pathway_output_ttest)<- NULL
  pathway_output <- join(pathway_output,pathway_output_ttest,type="full")
}

#Check p values within cell type for a concrete network. 
source_cells <- unique(reduced_files$source)
target_cells <- unique(reduced_files$target)
cell_cell_output<- data.frame()
for (i in source_cells){
  for(j in target_cells){
    try1 <- reduced_files %>% filter(source==i & target==j) %>% select(3:26)
    rownames(try1) <- try1$pathway_name
    try1$pathway_name <- NULL
    transpose_try1 <- t(try1)
    transpose_try1 <- data.frame(transpose_try1)
    Condition <- c("ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","ASD","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control")
    try2<-cbind(Condition,transpose_try1)
    rownames(try2) <- NULL
    try2<-try2[, which(colMeans(!is.na(try2)) > 0.3)]
    if (isTRUE(ncol(try2) >= 2)) {
      cell_cell_output_ttest <- list()
      pathways_ttest<-colnames(try2)
      pathways_ttest <- pathways_ttest[-1]
      for(x in pathways_ttest){
        fmla <- formula(paste0(x, " ~ Condition"))
        res<-wilcox.test(fmla, data = try2, exact = FALSE)
        cell_cell_output_ttest[[x]] <- res$p.value
      }
      cell_cell_output_ttest <- unlist(cell_cell_output_ttest)
      cell_cell_output_ttest<-data.frame(cell_cell_output_ttest)
      cell_cell_output_ttest$pathways<-rownames(cell_cell_output_ttest)
      rownames(cell_cell_output_ttest)<- NULL
      cell_cell_output_ttest <- cell_cell_output_ttest %>% filter(cell_cell_output_ttest < 0.05) %>% dplyr::rename(!!paste0(i,"_",j):=cell_cell_output_ttest)
      cell_cell_output <- join(cell_cell_output,cell_cell_output_ttest,type="full")
    }
  }
}

#Remove columns with all NA. At the end 141 cell-cell interactions with p <0.05 values
cell_cell_output2<-cell_cell_output[, colSums(is.na(cell_cell_output)) != nrow(cell_cell_output)]
#write.table(cell_cell_output2, "wilcoxon_pathway_by_cell_cell.xls", col.names = T, row.names = F, quote = F, sep="\t")

#Plot PCA of the samples using network probabilities.
#results_pathways[is.na(results_pathways)] = 0
rownames(summed_pathways_big) <- summed_pathways_big$pathways
summed_pathways_big$pathways <- NULL
summed_pathways_big_matrix <- data.matrix(summed_pathways_big)
t_summed_pathways_big_matrix<-t(summed_pathways_big_matrix)
t_summed_pathways_big_matrix<-t_summed_pathways_big_matrix[, which(colMeans(!is.na(t_summed_pathways_big_matrix)) > 0.3)]
t_summed_pathways_big_matrix[is.na(t_summed_pathways_big_matrix)] = 0
reduced_pca_pathways <- prcomp(t_summed_pathways_big_matrix, center = TRUE,scale. = TRUE)#not working with na
ggbiplot(reduced_pca_pathways,ellipse=TRUE,obs.scale = 1, var.scale = 1, var.axes = F, labels=rownames(t_summed_pathways_big_matrix), groups=Condition) +
  scale_colour_manual(name="Condition", values= c("#00C1AA", "red3"))+
  ggtitle("PCA based on pathways probability")+
  theme_minimal()+
  theme(legend.position = "bottom")

#Draw chord diagrams of some interesting networks. Calculate the mean values of ASD and Control probabilities for each network in each cell-cell pair and then do ASD-Control. In this case for NRG network.
test1 <- reduced_files %>% filter(pathway_name=="NRG") %>% select(!pathway_name)
test1$mean_asd <- apply(test1[,3:15],1,mean, na.rm=TRUE)
test1$mean_control <- apply(test1[,16:25],1,mean, na.rm=TRUE)
test1$mean_asd[is.na(test1$mean_asd)] = 0
test1$mean_control[is.na(test1$mean_control)] = 0
test1$mean_probability <- test1$mean_asd-test1$mean_control
test1 <- test1 %>% select (1,2,28)
#max(test1$mean_probability)
#min(test1$mean_probability)

#Since first I want to know which of the L-R are significant for chord diagram. For each interesting network
interaction_table <- readxl::read_xls("../all_interaction_wilcoxon_results.xls", na = "NA")
NRG_interactions <- interaction_table %>% filter(Pathways=="NRG") %>% select(!1&!2)
NRG_interactions<-NRG_interactions[, colSums(is.na(NRG_interactions)) != nrow(NRG_interactions)]
significant_NRG <- colnames(NRG_interactions)

#I need to create a data frame that consist in the soruce and targets and colors for each interaction
test1$concatenated <- str_c(test1$source,"_",test1$target)
#Once I have the significant interactions I will transform the values in 0.2 for positive change of signaling in asd and -0.2 for negative signaling in asd  in table test1. take into account that these values are arbitrary and when running chord diagram the strenght of the arrows will be based on that.
for (i in significant_NRG) {
  select1 <- test1 %>% filter(concatenated==i)
  if (isTRUE(select1$mean_probability > 0)){
    test1$concatenated[test1$concatenated %in% i] <-"#FF1515"
  }
}
for (i in significant_NRG) {
  select1 <- test1 %>% filter(concatenated==i)
  if (isTRUE(select1$mean_probability < 0)){
    test1$concatenated[test1$concatenated %in% i] <-"#476ECC"
  }
}
no_sig <- test1$concatenated
no_sig <- no_sig[!str_detect(no_sig,pattern="#")]
for (i in no_sig) {
  select1 <- test1 %>% filter(concatenated==i)
  if (isTRUE(select1$mean_probability > 0)){
    test1$concatenated[test1$concatenated %in% i] <-"#ffb2b2"
  }
}
for (i in no_sig) {
  select1 <- test1 %>% filter(concatenated==i)
  if (isTRUE(select1$mean_probability < 0)){
    test1$concatenated[test1$concatenated %in% i] <-"#C7D3EF"
  }
}

#Since I have only assigned until now the significant interactions:
#test1$concatenated[grep("_", test1$concatenated)] <- "#069943"
test2 <- test1 %>% select(1,2,3)
col_df <- test1 %>% select (1,2,4)
nameVals <- sort(unique(unlist(test2[1:2])))
myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
myMat[as.matrix(test2[c("source", "target")])] <- test2[["mean_probability"]]

#Draw chord diagram. Green means it is not significant. Red significant and up expresed in ASD and blue significant and down expressed in ASD
#Complete tutorial for chord diagram
#https://jokergoo.github.io/circlize_book/book/index.html

#If I would like to reorder the columns but keep colors.
#myMat<-symmpart(myMat)

#order <- c("AST-FB","Neu-NRGN-I","AST-PP","Endothelial","Neu-mat","IN-PV","IN-SST","IN-SV2C","IN-VIP", "L2/3","L4","L5/6","L5/6-CC","Oligodendrocytes","OPC","Neu-NRGN-II")


#myMat <- reorder_mat(mat=myMat, order = order)


#I can change names of the columns if I want to add a deg gene. And also I need to change them in the col_df dataframe.
#colnames(myMat)<-gsub("AST-FB","AST-FB\n(myMat)",colnames(myMat))
#rownames(myMat)<-gsub("AST-FB","AST-FB\n(myMat)",rownames(myMat))

#colnames(myMat)<-gsub("AST-PP","AST-PP\n(NRG2)",colnames(myMat))
#rownames(myMat)<-gsub("AST-PP","AST-PP\n(NRG2)",rownames(myMat))

#colnames(myMat)<-gsub("Endothelial","Endothelial\n(myMat)",colnames(myMat))
#rownames(myMat)<-gsub("Endothelial","Endothelial\n(myMat)",rownames(myMat))

#colnames(myMat)<-gsub("IN-PV","IN-PV\n(myMat)",colnames(myMat))
#rownames(myMat)<-gsub("IN-PV","IN-PV\n(myMat)",rownames(myMat))


#colnames(myMat)<-gsub("IN-SST","IN-SST\n(myMat)",colnames(myMat))
#rownames(myMat)<-gsub("IN-SST","IN-SST\n(myMat)",rownames(myMat))

#colnames(myMat)<-gsub("IN-SV2C","IN-SV2C\n(myMat)",colnames(myMat))
#rownames(myMat)<-gsub("IN-SV2C","IN-SV2C\n(myMat)",rownames(myMat))

#colnames(myMat)<-gsub("IN-VIP","IN-VIP\n(NRG1)",colnames(myMat))
#rownames(myMat)<-gsub("IN-VIP","IN-VIP\n(NRG1)",rownames(myMat))


#colnames(myMat)<-gsub("L2/3","L2/3\n(ERRB4)",colnames(myMat))
#rownames(myMat)<-gsub("L2/3","L2/3\n(ERRB4)",rownames(myMat))

#colnames(myMat)<-gsub("L4","L4\n(NRG2)",colnames(myMat))
#rownames(myMat)<-gsub("L4","L4\n(NRG2)",rownames(myMat))


#colnames(myMat)<-gsub("L5/6-CC","L5/6-CC\n(PTPRF)",colnames(myMat))
#rownames(myMat)<-gsub("L5/6-CC","L5/6-CC\n(PTPRF)",rownames(myMat))

#colnames(myMat)[10]<-gsub("L5/6","L5/6\n(myMat)",colnames(myMat)[10])
#rownames(myMat)[10]<-gsub("L5/6","L5/6\n(myMat)",rownames(myMat)[10])

#colnames(myMat)<-gsub("Microglia","Microglia\n(myMatAP2)",colnames(myMat))
#rownames(myMat)<-gsub("Microglia","Microglia\n(myMatAP2)",rownames(myMat))

#colnames(myMat)<-gsub("Neu-mat","Neu-mat\n(LRRC4C)",colnames(myMat))
#rownames(myMat)<-gsub("Neu-mat","Neu-mat\n(LRRC4C)",rownames(myMat))

#colnames(myMat)[14]<-gsub("Neu-NRGN-I","Neu-NRGN-I\n(myMat)",colnames(myMat)[14])
#rownames(myMat)[14]<-gsub("Neu-NRGN-I","Neu-NRGN-I\n(myMat)",rownames(myMat)[14])

#If I would like to change the order I would need to add the new colum number
#colnames(myMat)[15]<-gsub("Neu-NRGN-II","Neu-NRGN-II\n(myMat)",colnames(myMat)[15])
#rownames(myMat)[15]<-gsub("Neu-NRGN-II","Neu-NRGN-II\n(myMat)",rownames(myMat)[15])

#colnames(myMat)<-gsub("Oligodendrocytes","Oligo",colnames(myMat))
#rownames(myMat)<-gsub("Oligodendrocytes","Oligo",rownames(myMat))

#colnames(myMat)<-gsub("OPC","OPC\n(NRGAP2)",colnames(myMat))
#rownames(myMat)<-gsub("OPC","OPC\n(NRGAP2)",rownames(myMat))

#col_df$source[col_df$source %in% "L4"] <-"L4\n(NRG2)"
#col_df$target[col_df$target %in% "L4"] <-"L4\n(NRG2)"

#col_df$source[col_df$source %in% "L2/3"] <-"L2/3\n(ERRB4)"
#col_df$target[col_df$target %in% "L2/3"] <-"L2/3\n(ERRB4)"

#col_df$source[col_df$source %in% "AST-PP"] <-"AST-PP\n(NRG2)"
#col_df$target[col_df$target %in% "AST-PP"] <-"AST-PP\n(NRG2)"

#col_df$source[col_df$source %in% "Oligodendrocytes"] <-"Oligo"
#col_df$target[col_df$target %in% "Oligodendrocytes"] <-"Oligo"

#Initialize grid colors
grid.col = NULL
order <- c("AST-FB","Neu-NRGN-I","AST-PP","Endothelial","Neu-mat","IN-PV","IN-SST","IN-SV2C","IN-VIP", "L2/3","L4","L5/6","L5/6-CC","Oligodendrocytes","OPC","Neu-NRGN-II")
# Set row and column labels colors
grid.col[colnames(myMat)] = c("#E41A1C", "#377EB8",  "#984EA3", "#4DAF4A","#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", "#E7298A","#910241", "#00CDD1")
grid.col[rownames(myMat)] = c("#E41A1C", "#377EB8",  "#984EA3", "#4DAF4A","#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", "#E7298A","#910241", "#00CDD1")
# Parameters for circos layout. The gap.degree specifies the gap between two neighbour sectors. It can be a simyMate value or a vector. If it is a vector, the first value corresponds to the gap after the first sector
circos.par(gap.degree = 4)

#I have a matrix but it can be woorked too with data frame:
#col_mat = myMat
#dim(col_mat) = dim(myMat)#VIP
#Sometimes depedning on the strenght values I will have to change these values
#Now we can assign a concrete color to the values we want as follow:
#col_mat[myMat == 0.2] = "#FF1515"
#col_mat[myMat == -0.2] = "#476ECC"
#col_mat[myMat < 0.2 & myMat > -0.2] = "#069943"


#It can only be used with data frame


#col_mat = function(x) ifelse(x < 0.2 & x > -0.2, "#069943", ifelse(x==0.2, "#FF1515","#476ECC"))

#If I add this one no longer lines same color as grids.
chordDiagram(myMat, grid.col = grid.col, col = col_df,directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", annotationTrack = "grid", scale=F,preAllocateTracks = list(list(track.height = 0.01), list(track.height = 0.01)))

#annotationTrackHeight = mm_h(5)
#I use annotationTrackHeight = mm_h(5) #using this one I can add the names of the cells in the grid.          


#Scale can be turned off or not.


#To add names to the grid.
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), mean(ylim), sector.index, facing = "inside", niceFacing = TRUE, cex=1)
}, bg.border = NA)

#For color: col="white"

#Other options for facing: "bending.inside", "inside"
title("NRG network", cex = 1)
circos.clear()

