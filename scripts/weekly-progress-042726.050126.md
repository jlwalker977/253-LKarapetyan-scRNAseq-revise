Progress for the week of 4.27.26
4.27
TIL-redo-recleaned obj completed and was ready to use
Sent the metadata (from TIL only) to Arjun for his ref
Tested the percentage of IG genes present in the TIL and PBMC populations at each cleaning step
FeaturePlots of genes representative of Macs, CD3, CD8, CD4 T cells, and 4x IG genes that were detected in previous deg runs
Built slide deck of these
Sent to Brian to address his concerns from a previous email

4.28
Quick slingshot trajectories for the newly cleaned objects, separated into their general cell types and downsampled 500cells/cluster
The trajectories will be on the cell_types2 (the more specific ones from Arjun originally), so that they're informative
Ran Sling on PCA and viewed it in PCA
Ran Sling on UMAP and viewed it in UMAP
Found that running it in PCA and mapping the coords onto UMAP (recommended in teh literature), was overly complicated for our purposes. If Lilit wants to do more of it, I will do it in a more precise way
Collected refs to read up on:
IG genes expressed in T cells as non-coding transcripts
Can be missing start codon so they never get made into proteins

4.29
Got email from Adam Soloff's group that the issue with their rnaseq files was corrected (I saw that there was one folder that was way smaller than the others)
Started the new Globus download and let it run overnight
Wrote Globus SOP for sending and receiving
253 - DEG loops on cleaning objects, PBMC and TIL
Ran them as background jobs
Some debugging with Claude led to better heatmap annotation

4.30
Combined the DEG plots from 4.29 into a powerpoint
Volcanoes on top, heatmaps on the bottom
Ran recurring-deg-barplot functions with Claude, marked with directionality
Combined TIL and PBMC metadata csvs to send to Arjun

5.1
Coursera work on IBM R course (due to finish it by June and only 2.5/5 right now)
Background reading for projects 253 and 260
Delete project 260 folders off of /ix/ since they're all moved over to /vast/ now
