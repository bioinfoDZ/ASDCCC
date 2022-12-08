"Script 1: Global approach CellChat"
#CellChat analysis summary#

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(data.table)
library(tidyverse)
library(plyr)
library(dplyr)
library(NMF)
library(ComplexHeatmap)
library(ggalluvial)
library(vioplot)


#GLOBAL ANALYSIS IN WHICH ALL COMPUTED VALUES ARE BASED IN THE COMBINATION OF SAMPLES WITHIN SAME CONDITION. 
#Matrix from the paper was split in Control and ASD samples. 

#Create CellChat object using metadata and data matrix. Colnames in data matrix have to be the same as rownames in metadata.
Meta_data <- read.table("/Users/mastorkia/Dropbox (EinsteinMed)/scRNASeq_Projects/AUTISM/Autism_DATA/meta.tsv", sep="\t", header=TRUE) 
ASD_pre_cortex_meta_data<-Meta_data%>%filter(diagnosis=="ASD")%>%filter(region=="PFC")
ASD_pre_cortex_meta_data<- ASD_pre_cortex_meta_data%>% column_to_rownames(var="cell")

#This matrix was split by the previous main matrix with all control and ASD patients. 
ASD_matrix <- fread("/Users/mastorkia/Dropbox (EinsteinMed)/scRNASeq_Projects/AUTISM/Autism_DATA/ASD_names_matrix.csv")
ASD_matrix<-ASD_matrix %>% distinct(external_gene_name, .keep_all = TRUE)
ASD_matrix <- ASD_matrix %>% column_to_rownames(var="external_gene_name")

#Filter the matrix too. 
keep_prefrontal_cortex_cells<-rownames(ASD_pre_cortex_meta_data)
ASD_pre_cortex_matrix <- ASD_matrix[, keep_prefrontal_cortex_cells]
ASD_pre_cortex_matrix <- data.matrix(ASD_pre_cortex_matrix)

#Do the same with Control samples.


#Create the CellChat object for ASD and Control separately.
ASD_pre_cortex_cellchat <- createCellChat(object = ASD_pre_cortex_matrix, meta = ASD_pre_cortex_meta_data, group.by = "cluster")
#Get number of cells for each cell type:
groupSize_asd_pre_cortex <- as.numeric(table(ASD_pre_cortex_cellchat@idents))
idents_asd_pre_cortex<-levels(ASD_pre_cortex_cellchat@idents)
asd_pre_cell_type_number<-cbind(groupSize_asd_pre_cortex,idents_asd_pre_cortex)
CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
ASD_pre_cortex_cellchat@DB <- CellChatDB.use
ASD_pre_cortex_cellchat <- subsetData(ASD_pre_cortex_cellchat) # subset the expression data of signaling ligands for saving computation cost
ASD_pre_cortex_cellchat <- identifyOverExpressedligands(ASD_pre_cortex_cellchat)
ASD_pre_cortex_cellchat <- identifyOverExpressedInteractions(ASD_pre_cortex_cellchat)
ASD_pre_cortex_cellchat <- projectData(ASD_pre_cortex_cellchat, PPI.human)

#Compute interaction probabilities
ASD_pre_cortex_cellchat <- computeCommunProb(ASD_pre_cortex_cellchat, raw.use = TRUE, type = "truncatedMean", trim = 0.2)#The average gene expression is zero if the percent of expressed cells in one group is less than 20%
ASD_pre_cortex_cellchat <- filterCommunication(ASD_pre_cortex_cellchat, min.cells = 10)
ASD_pre_cortex_cellchat <- computeCommunProbPathway(ASD_pre_cortex_cellchat)
ASD_pre_cortex_cellchat <- aggregateNet(ASD_pre_cortex_cellchat)

#Same for Control object. 

#From the created objects get the probabilities of Ligand-Receptors (L-R) (net) or networks (netP) in which L-R are involved and save those to use it later. 
ASD_pre_net <- subsetCommunication(ASD_pre_cortex_cellchat)
#write.table(ASD_pre_net,"asd_pre_cortex_net.csv", quote=F, sep="\t", row.names = F)
Control_pre_net <- subsetCommunication(Control_pre_cortex_cellchat)
#write.table(Control_pre_net,"control_pre_cortex_net.csv", quote=F, sep="\t", row.names = F)

ASD_pre_netP <- subsetCommunication(ASD_pre_cortex_cellchat, slot.name = "netP")
#write.table(ASD_pre_netP,"asd_pre_cortex_netP.csv", quote=F, sep="\t", row.names = F)
Control_pre_netP <- subsetCommunication(Control_pre_cortex_cellchat, slot.name = "netP")
#write.table(Control_pre_netP,"control_pre_cortex_netP.csv", quote=F, sep="\t", row.names = F)

#This same analysis willl be repeated for each sample separately so we can run the "Sample-by-sample approach" with the information of L-R and networks. 

#Continuing with "Global approach"
#Draw chord diagrams of cell-cell interactions in Control and ASD datasets
groupSize_Control <- as.numeric(table(Control_pre_cortex_cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(Control_pre_cortex_cellchat@net$count, vertex.weight = groupSize_Control, weight.scale = T, label.edge= F, title.name = "Number of interactions in Control dataset")
netVisual_circle(Control_pre_cortex_cellchat@net$weight, vertex.weight = groupSize_Control, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


#Check outgoing and incoming network signaling within each dataset.
Control_pre_cortex_cellchat <- netAnalysis_computeCentrality(Control_pre_cortex_cellchat, slot.name = "netP") 
ASD_pre_cortex_cellchat <- netAnalysis_computeCentrality(ASD_pre_cortex_cellchat, slot.name = "netP") 

#Identify signaling patterns in order to see how cells talk together.First run two metrics (Cophenetic and Silhouette) that are based in hierachical clustering. To check the plot just type k2_control.
k2_control<-selectK(
  Control_pre_cortex_cellchat, #ASD_pre_cortex_cellchat
  slot.name = "netP",
  pattern = "incoming", #"outgoing"
  title.name = "Infer Incoming Control patterns",
  do.facet = TRUE,
  k.range = seq(2, 10),
  nrun = 30,
  seed.use = 10
)

#Then based on the plot set the number of patterns for each outgoing, incoming Control, ASD.  
Control_pre_cortex_cellchat <- identifyCommunicationPatterns(
  Control_pre_cortex_cellchat,
  slot.name = "netP",
  pattern = "outgoing",
  k = 5
)

#Plot the patterns in a river plot
river_control <- netAnalysis_river(
  Control_pre_cortex_cellchat,
  cutoff = 0.5, #they use 0.5 threshold in order to plot the most enriched patterns and networks
  slot.name = "netP",
  pattern = "outgoing",
  color.use.signaling = "grey50",
  do.order = FALSE,
  main.title = "Outgoing signaling Control",
  font.size = 2.5,
  font.size.title = 12
)

#Now lets plot with dots. Data for this plot:
dots_control <- netAnalysis_dot(
  Control_pre_cortex_cellchat,
  slot.name = "netP",
  pattern = "outgoing",
  shape = 21,
  dot.size = c(1, 5),
  dot.alpha = 1,
  main.title = "Outgoing signaling patterns-Control",
  font.size = 10,
  font.size.title = 12
)



#Get the information of the patterns and save it.
Control_outgoing_cell_pattern <- data.frame(Control_pre_cortex_cellchat@netP$pattern[["outgoing"]][["pattern"]][["cell"]])
Control_outgoing_siganling_pattern <- data.frame(Control_pre_cortex_cellchat@netP$pattern[["outgoing"]][["pattern"]][["signaling"]])
Control_incoming_siganling_pattern <- data.frame(Control_pre_cortex_cellchat@netP$pattern[["incoming"]][["pattern"]][["signaling"]])
Control_incoming_cell_pattern <- data.frame(Control_pre_cortex_cellchat@netP$pattern[["incoming"]][["pattern"]][["cell"]])

#write.table(Control_outgoing_cell_pattern, "Control_outgoing_cell_pattern.csv", quote = F, sep="\t", row.names = F)
#write.table(Control_incoming_cell_pattern, "Control_incoming_cell_pattern.csv", quote = F, sep="\t", row.names = F)
#write.table(Control_outgoing_siganling_pattern, "Control_outgoing_signaling_pattern.csv", quote = F, sep="\t", row.names = F)
#write.table(Control_incoming_siganling_pattern, "Control_incoming_signaling_pattern.csv", quote = F, sep="\t", row.names = F)

ASD_outgoing_cell_pattern <- data.frame(ASD_pre_cortex_cellchat@netP$pattern[["outgoing"]][["pattern"]][["cell"]])
ASD_outgoing_siganling_pattern <- data.frame(ASD_pre_cortex_cellchat@netP$pattern[["outgoing"]][["pattern"]][["signaling"]])
ASD_incoming_siganling_pattern <- data.frame(ASD_pre_cortex_cellchat@netP$pattern[["incoming"]][["pattern"]][["signaling"]])
ASD_incoming_cell_pattern <- data.frame(ASD_pre_cortex_cellchat@netP$pattern[["incoming"]][["pattern"]][["cell"]])

#write.table(ASD_outgoing_cell_pattern, "ASD_outgoing_cell_pattern.csv", quote = F, sep="\t", row.names = F)
#write.table(ASD_incoming_cell_pattern, "ASD_incoming_cell_pattern.csv", quote = F, sep="\t", row.names = F)
#write.table(ASD_outgoing_siganling_pattern, "ASD_outgoing_signaling_pattern.csv", quote = F, sep="\t", row.names = F)
#write.table(ASD_incoming_siganling_pattern, "ASD_incoming_signaling_pattern.csv", quote = F, sep="\t", row.names = F)


#Until now everything has been using each object separately. Now merge both and work with a unique object.
object.list <- list(Control = Control_pre_cortex_cellchat, ASD = ASD_pre_cortex_cellchat)
pre_cortex_cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#Group one will be Control, group 2 ASD
gg1 <- compareInteractions(pre_cortex_cellchat, show.legend = F, group = c(1,2), width = 0.4)
gg2 <- compareInteractions(pre_cortex_cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(pre_cortex_cellchat, weight.scale = T, title.name = "Differential number of interaction", arrow.size=0.1, arrow.width=0.5)
netVisual_diffInteraction(pre_cortex_cellchat, weight.scale = T, measure = "weight", arrow.size=0.1, arrow.width=0.5)

#In order to plot counts
ASD_counts<-ASD_pre_cortex_cellchat@net$count
Control_counts<-Control_pre_cortex_cellchat@net$count
diff_counts<-ASD_counts-Control_counts

#All L-R  edge strengths
ASD_weight<-ASD_pre_cortex_cellchat@net$weight
Control_weight<-Control_pre_cortex_cellchat@net$weight
diff_weight<-ASD_weight-Control_weight

#Mixed heatmap. In order to see better the changes and also to have the row with same order. 
col_fun = circlize::colorRamp2(c(-30,0, 30), c("blue", "gray", "red"))
Heatmap(diff_counts, name="ASD vs Control pre Cortex Counts", rect_gp = gpar(col = "white", lwd = 2), column_dend_reorder=TRUE, column_title = "Counts", cluster_columns = FALSE, col = col_fun) + 
Heatmap(diff_weight, name="ASD vs Control pre Cortex Strengh", rect_gp = gpar(col = "white", lwd = 2), column_dend_reorder=TRUE, column_title = "Strengh", cluster_columns = FALSE)



#Plot the probability values
col_weight = circlize::colorRamp2(c(-0.5,0, 0.5), c("blue", "gray", "red"))
strengh_heatmap <- Heatmap(diff_weight, name="Strengh", rect_gp = gpar(col = "white", lwd = 2), column_dend_reorder=TRUE, cluster_columns = FALSE, cluster_rows = FALSE,col = col_weight, row_title="Sending cells", column_title="Receiving cells", column_title_side = "bottom",
                           column_title_gp = gpar(fontsize = 14),row_title_side =  "right",
                           row_title_gp = gpar(fontsize = 14), cell_fun = function(j, i, x, y, width, height, fill) {
                             grid.text(sprintf("%.2f", diff_weight[i, j]), x, y, gp = gpar(fontsize = 6))
                           }
)
strengh_heatmap


#Draw chord diagram with the number of interactions in ASD and Control
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#Identify differentially expressed L-R on networks based on the sum of all cell-cell communication probabilities for that L-R or sum of the L-R involved in network. 
#Here it uses the Wilcoxon non parametrical test to test between two conditions.
#From here we can extract the down and up regulated networks for ASD. As well as down-up L-R interactions. 
Significant_differences <- rankNet(
  pre_cortex_cellchat,
  slot.name = "netP", #If I type net it will do it using L-R relations
  mode = "comparison",
  comparison = c(1, 2),
  stacked = T,
  do.stat = T,
  cutoff.pvalue = 0.05,
  tol = 0, #tolerance can be key. so decide where to cut it
  thresh = 0.05,
  show.raw = T,
  return.data = T,
  font.size = 9,
  title = "Network differences between ASD-Control conditions",
  bar.w = 0.60, #the width of bar plot
  rel_heights = c(0.9, 0, 0.1),
  ylim = c(0, 50)
)

#write.table(Significant_differences$signaling.contribution, "significant_networks_pre_cortex.csv", quote=F, sep="\t", row.names = F)```


