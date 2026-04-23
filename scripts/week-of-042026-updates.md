04/23/26
Daily time card
Started at 1030, overslept, will stay on until 630 to make up the time
Waiting to hear back from Brian about the scDblFinder run and if I shold keep it in the cleaned object or not
Asked more specifics about the recurring DEG table and how best to represent it
Asked for his reasoning on things
Reading ASoloff's new project info he just sent over
Brushing up on Immunology background info so I can better understand the results I get

4/22/26
Daily time card
Checked output from DEG run on clean4-TIL-object 
Saw IG genes in CD3 and CD4 cell types where they are not supposed to be
Looked into how to better clean the B cell contamination
Chatted with Claude about DoubletFinder within Seurat
Decided to go with scDblFinder instead bc it had better Seurat integration
Ran scDblFinder
It requires separating objects by sample-id, for us it was Patient.Number
It reran PCA and UMAP so the UMAP looks different now
There is one little CD4 T cell protrusion on the macrophage cluster, how to remove it?
Tried to determine what fine-cell-type it was
Checked a bunch of gene expression using featurePlot for the subtypes: Early activated/memory, gamma delta T cells, and ISG?macs
Inconclusive
Put all of this in an email to Brian to ask what to proceed with and what to shelve
Slingshot runs from 4/21 were unsuccessful, trying to determine reason
Two were noted as failing, the other 4 never showed any output
Will need to put more diagnostics in to figure out why

4/21/26
Daily time card
Attempted SoupX and DoubletFinder runs, unsuccessful
Brian said to just remove the clump of cells manually
removed clusters 6+9
Still saw gene expression patterns that were troubling
Ran decontX instead of SoupX bc its better with aggregated seurat objects and doesn't require the raw data
Ran cleaned object through DEG loop (TIL) to check to see if the IG genes were still there
Started slingshot runs on indiv cell types (looking at fine-cell-types within those general ones) for PBMC and TIL as background jobs

4/20/26
Daily time card
Recurring DEG directionality
Talked to Claude about what to do to show this in a visual way (preferrably heatmap)
Created heatmaps (annotated) of the direction of recur degs separated by comparison type
Asked claude how to interpret
Claude determined that the TIL PFS and timepoint plots showed good TIL signal
PBMC plots were less strong
Emailed them to Brian and he didn't understand how/why I had made them, tried to explain reasoning
Ran TIL and PB through DEG loop to make sure the gene filtering is working