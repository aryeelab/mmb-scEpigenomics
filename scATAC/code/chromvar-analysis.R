# Credit to Alicia Schep 
# for developing a v1 of this script

library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19) # change based on reference genome
library(SummarizedExperiment)
library(chromVARmotifs)

# Initialize parallel processing
library(BiocParallel)
register(MulticoreParam(2)) # adjust according to your machine

# Import/filter data (replace with appropriate file paths)
# Pre-processing already finished
if(FALSE){
  peakfile <- "data/peaks.bed"
  peaks <- getPeaks(peakfile)
  bamfiles <- list.files("data/bams", full.names = TRUE)
  raw_counts <- getCounts(bamfiles, peaks, paired =  TRUE, by_rg = FALSE,
                          format = "bam", colData = DataFrame(source = bamfiles))
  
  # Filter low quality samples and peaks
  counts_filtered <- filterSamples(raw_counts, min_depth = 500,
                                   min_in_peaks = 0.15, shiny = FALSE)
  counts <- filterPeaks(counts_filtered)
}
load("../data/example_counts.rda")
counts <- example_counts

# Get GC content/peak; get motifs from chromVARmotifs package; find kmers
counts <- addGCBias(counts, genome = BSgenome.Hsapiens.UCSC.hg19)
data("human_pwms_v2") # also mouse_pwms_v2
motif_ix <- matchMotifs(human_pwms_v2, counts, genome = BSgenome.Hsapiens.UCSC.hg19)
dim(motif_ix)

# Compute deviations; typically most time consuming step
dev <- computeDeviations(object = counts, annotations = motif_ix)

# Find variable motifs
variabilityAll <- computeVariability(dev)
plotVariability(variabilityAll, use_plotly = FALSE)

# Visualize single cells with a tSNE
tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 10,
                               shiny = FALSE)
tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation = "CTCF", shiny = FALSE)
tsne_plots