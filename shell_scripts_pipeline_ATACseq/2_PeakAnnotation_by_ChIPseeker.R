# This script performs genomic annotation on selected peaks using ChIPseeker. Its output is input to clusterProfiler.
# Input file can be of any format but must have 3 columns of bed format. Chromosomes must have chr because of
# TxDb.Hsapiens.UCSC.hg38.knownGene

require(dplyr)
require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(EnsDb.Hsapiens.v86)
require(AnnotationDbi)

# STEP 1: Loading data ####

smpls_ = as.list(list.files("../out", pattern= "tsv", full.names=T))

names(smpls_) = lapply(X = smpls_, FUN = function(fl_) { return(tail(strsplit(x = fl_, split = "[/]")[[1]],1)) })

# STEP 2: Annotating ####

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene      # assigning UCSC annotation database

peakAnnoList = lapply(smpls_, FUN = function(smpl_){ annotatePeak(peak = smpl_, TxDb = txdb, tssRegion = c(-1000, 1000), verbose = T) })

for(fl_ in names(smpls_))
{
  annot_ = data.frame(peakAnnoList[[fl_]]@anno)
  
  # matching ENTREZ ID to HGNC gene symbol
  entrez = annot_$geneId                                                                                                        # Get the entrez IDs
  annotations_edb = AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = entrez, columns = c("GENENAME"), keytype = "ENTREZID")     # Return the gene symbol for the set of Entrez IDs
  annotations_edb$ENTREZID = as.character(annotations_edb$ENTREZID)                                                             # Return the gene symbol for the set of Entrez IDs
  annot_ = left_join(x = annot_, y = annotations_edb, by = c("geneId"="ENTREZID"))
  
  write.table(x = annot_, file = paste0("../out/", strsplit(x = fl_, split = '[.]')[[1]][1] ,"_annot.tsv"), sep="\t", quote = F, row.names = F)
}

# STEP 3: Plotting annotation statistics ####

pdf(file = '../out/annotations_plots.pdf', width = 10)
plotAnnoBar(peakAnnoList)                                                                 # bar chart of genomic feature representation
plotDistToTSS(peakAnnoList, title="Distribution of binding loci \n relative to TSS")      # distribution of binding site relative to TSS
graphics.off()
