#Outputs from ASD_Cellchat_Script3#

As in the script 2, here I call all the files with all L-R interactions in each sample.

As before I create a big table with all L-R interaction probabilities for each sample.

I follow running Wilcoxon rank test to identify specific L-R interactions between cell-cell.

Also I calculate the fold change of those interactions therefore I know if the significant interactions are up or down regulated in our conditions. 

I code for plotting these results in a big Heatmap. In which I show in y axes the L-R interactions and in X axes the cell-cell interactions. Color determine up or down. 

Using all p values (not filtered) I plot a value histogram.

I also compute FDR multiple testing.

And finally I trace back the pathway name of those significant interactions. 