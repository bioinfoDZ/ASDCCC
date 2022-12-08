#Outputs from ASD_Cellchat_Script1#

#Global analysis#

First steps describe how to load and manage the data matrix and ork with metadata and matrix to create the CellChat object.

Then following CellChat vignette create the CellChat object for ASD and Control samples. In this case all samples from ASD will be analyzed together, same for Control samples.

Using subsetCommunication (slot.name="netP" or "net") get all the interactions for Control and ASD as well as pathways interactions and save it as excel files.

This same analysis can be done for each sample. Net and netP results for each sample will be saved as well for further "Sample-by-sample analysis".

With the object and following the vignette get the main figures:

Chord diagrams of cell-cell interactions.

Centrality analysis to get outgoing and incoming signaling as well as clustering of cells based on sending pathways or receiving pathways. 

After merging Control and ASD objects:

Heatmaps with difference in counts and strength between both conditions.

Chord diagrams reflecting up and down interactions.

And lastly significantly different pathways between the two conditions as well as relative contribution of the pathways in each of the samples.