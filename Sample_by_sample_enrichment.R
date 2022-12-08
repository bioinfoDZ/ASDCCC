"Script4: Sample by sample enrichment analysis"

#With the genes involved in significant L-R interactions we do an enrichment analysis.
library(org.Hs.eg.db) 
library(DOSE)
library(enrichplot)
library(ggnewscale)
library(clusterProfiler)
library(ggplot2)
library(GeneOverlap)
library(biomaRt) #2.44.4
library(openxlsx)
library(dplyr)
library(plyr)


#Run GO enrichment
#These are the genes that are involed in the significant interactions from script 3. 
all_diff_network_genes <- read.table("../enrichment_genes.txt")
all_diff_network_genes<-all_diff_network_genes$V1
ego <- enrichGO(all_diff_network_genes,
                OrgDb = "org.Hs.eg.db", 
                ont="BP",
                keyType = "SYMBOL", 
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01
)
#remove redundant pathways  
ego2 <- simplify(ego)
go_bp_all_genes<- data.frame(ego2@result)


#Plot the first 30 significant networks.
reduced_go_bp_all_genes<- go_bp_all_genes[1:30,]
reduced_go_bp_all_genes$'-log10(FDR.p.value)'<- -log10(reduced_go_bp_all_genes$p.adjust)
ggplot( reduced_go_bp_all_genes, aes( x=0.05, y=Description)) +  
  geom_dotplot( binaxis = "y", stackdir = "center", binwidth=0.1)+ 
  scale_x_continuous(NULL, breaks = NULL) + coord_cartesian(xlim=c(0, 0.1)) + 
  geom_point(aes(size=Count, colour=`-log10(FDR.p.value)`)) + ylab (NULL)

#Now run an enrichment analysis based on knwon genes involved in close diseases. For that run Fisher test.
#According to GENCODE's last release there are 19951 protein coding genes.
protein_coding_genes <- 19951

#Using a table created with gene sets related to neurodevelopmental disorders...
diseases_gene_sets <- read.xlsx("../disease_gene_sets.xlsx")
tested_overlap_results <- list()
for (x in colnames(diseases_gene_sets)){
  overlap_results <- newGeneOverlap(all_diff_network_genes, 
                                    diseases_gene_sets[[x]], 
                                    protein_coding_genes)
  tested_overlap_results[[x]]<- testGeneOverlap(overlap_results)
}

#Plot results from diseases enrichment analysis. 
enrichment_results <- read.xlsx("../enrichemnt_disease_results.xlsx")
enrichment_results$'-log10(p value)' <- -log10(enrichment_results$p.value)
#reorder by logp
enrichment_results$Name <- reorder(enrichment_results$Name, enrichment_results$`-log10(p value)`)

h <- ggplot(enrichment_results, aes( x=0.05, y=Name)) +  
  geom_dotplot( binaxis = "y", stackdir = "center", binwidth=0.1)+ 
  scale_x_continuous(NULL, breaks = NULL) + coord_cartesian(xlim=c(0, 0.1))+
  geom_point(aes(size=Counts, colour=`-log10(p value)`)) + ylab (NULL)

h

