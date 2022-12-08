#Outputs from ASD_Cellchat_Script2#

#Start with Sample by sample approach#

This is the first step to run sample-by-sample steps. 

We already have run CellChat for each sample and obtained files for all the L-R interactions and for all pathways interactions for each sample.

In this script we will call the files with pathway (netP) information and merge them in a big table with interaction strengths of all pathways and samples. 

And from here we will run Wilcoxon rank test in 2 ways:

1-Taking into account pathway probability in all the cell-cell interactions together. So summing all cell-cell interactions for each pathway.

2- Taking into account pathway probability for each cell-cell interaction. 

Also we will run principal components based on pathway probabilities.

Last code is to draw chord diagram for each of the pathways in which difference between conditions is portrait with different colors and arrows show the interaction directions (up or down). 