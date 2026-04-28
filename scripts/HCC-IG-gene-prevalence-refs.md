This is an area without a single firmly established benchmark, but there are reasonable reference points from the literature:

**Published ranges in tumor-infiltrating T cells**
- Most studies report somewhere between **5-15%** of T cells showing detectable IG gene expression in scRNAseq
- Studies specifically in HCC (which is your context) tend to report on the higher end due to the immune-rich TME and frequent tertiary lymphoid structure presence
- Your PBMC CD3plus at ~12% sits right in the normal range
- Your TIL CD3plus at ~23% is elevated but not implausible for HCC

**What drives variation across studies:**
- Sequencing depth — deeper sequencing captures more low-level transcripts, inflating detection rates
- Dissociation protocol — longer dissociation times increase ambient RNA release from plasma cells
- 10x chemistry version — v3 captures more transcripts per cell than v2, increasing ambient signal detection
- Whether ambient RNA correction was applied before reporting
- How the T cell cluster was defined — broad CD3plus vs stringent CD8/CD4 subclusters will differ substantially

**HCC-specific context**
- HCC is notable for having significant B cell and plasma cell infiltration, often organized into tertiary lymphoid structures
- Several HCC scRNAseq papers (Zheng et al. 2017 Cell, Zhang et al. 2019 Cell) report detectable IG transcripts in T cell clusters and attribute them to a mix of ambient RNA and genuine low-level transcription
- The consensus in those papers is that IG+ T cells below ~25% detection with no accompanying constant region expression are not considered a quality concern

**The more informative comparison than a fixed threshold** is your own internal controls:
- Your Macrophages at ~12% TIL serve as a useful baseline — macrophages have no biological reason to express IG genes, so their rate approximates your ambient RNA floor
- Your CD3plus at ~23% TIL is roughly 2x that floor, which suggests a mix of ambient and genuine transcription
- Your PBMC CD3plus at ~12% essentially matches the macrophage ambient floor, suggesting PBMC T cells have minimal genuine IG transcription as expected for circulating cells

**Bottom line for your data:**
The TIL elevation is within a biologically plausible range for HCC specifically, the cleaning steps didn't meaningfully change it, and constant region co-expression would be the most important next check. If `IGKC`/`IGLC`/`IGHM` are absent in those IG-positive CD3plus cells, the ~23% figure is defensible as non-productive transcription plus ambient background rather than a quality problem.

Google scholar search terms:
# HCC scRNAseq T cell papers
"hepatocellular carcinoma" "single cell" "tumor infiltrating lymphocytes" 
"T cell" "scRNA-seq"

# Ambient RNA and IG genes specifically  
"ambient RNA" "scRNA-seq" "immunoglobulin" "T cells"
"decontX" OR "SoupX" "immunoglobulin contamination"

# Non-productive rearrangement in T cells
"non-productive rearrangement" "T cells" "immunoglobulin transcripts"
"sterile transcription" "immunoglobulin" "T lymphocytes"

# General TIL scRNAseq quality
"tumor infiltrating lymphocytes" "single cell" "quality control" 
"ambient RNA" "doublets"