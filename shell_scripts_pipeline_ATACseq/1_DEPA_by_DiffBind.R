# This script performs peak/binding site enrichment analysis using DiffBind. Its output is input to ChIPseeker

require(DESeq2)
require(edgeR)
require(DiffBind)

alpha_ = 0.01     # significance level

# STEP1: Performing DEPA ####

conds_ = c('IFNy','vehicle','CD28-negative','CD28-positive')
for(cont_ in combn(length(conds_),2,simplify = F) )
{
  fl_ = paste(conds_[cont_],collapse = '_vs_')
  cat('Processing ',fl_,'\n')
  
  smpls_ = tryCatch(expr = suppressWarnings(read.csv(file = paste0('../data/samplesheet_',fl_,'.csv'), stringsAsFactors = T)), error = function(er) {return(NULL)})
  if(is.null(smpls_)){ next() }
  
  isFIle = tryCatch(expr = suppressWarnings(load(file = paste0('../data/DiffBindObj_',fl_,'.RData'))), error = function(er) {return(NULL)})
  if(is.null(isFIle))
  {
    dba_ = dba(sampleSheet = smpls_)                                          # constructing the consensus peak set (union set) covering all peak regions in all conditions
    dba_ = dba.count(dba_, bUseSummarizeOverlaps=TRUE)                        # constructing binding affinity matrix (read count matrix)
    dba_ = dba.contrast(dba_, categories = DBA_CONDITION, minMembers = 2)     # creating contrasts for experiments w/ replicates
    dba_ = dba.analyze(dba_, method = DBA_ALL_METHODS)                        # differential peak enrichment analysis (or differential binding site affinity/enrichment analysis)
    save(dba_, file = paste0('../data/DiffBindObj_',fl_,'.RData'))
  }
  
  # STEP2: Writing DEPA result ####
  message('plotting DPEA result')
  
  # plotting DPEA result
  pdf(file = paste0('../out/',fl_,'.pdf'), width = 8, height = 8)
  
  tryCatch(expr = dba.plotPCA(dba_, contrast = 1, th = alpha_, bUsePval = F, method = DBA_EDGER, attributes = DBA_CONDITION, label = DBA_ID, sub = 'edgeR'), error = function(e_) { return(NA) })      # PCA on significant DE peaks (reducing dimensions created by all DE peaks)
  tryCatch(expr = dba.plotPCA(dba_, contrast = 1, th = alpha_, bUsePval = F, method = DBA_DESEQ2, attributes = DBA_CONDITION, label = DBA_ID, sub = 'DESeq'), error = function(e_) { return(NA) })
  
  tryCatch(expr = dba.plotMA( dba_, contrast = 1, th = alpha_, bUsePval = F, method = DBA_EDGER, bLoess = F, fold = 2, sub = 'edgeR'), error = function(e_) { return(NA) })                                                  # MA plot for edgeR called DE peaks (the distribution of up/down-regulated peaks)
  tryCatch(expr = dba.plotMA( dba_, contrast = 1, th = alpha_, bUsePval = F, method = DBA_DESEQ2,bLoess = F, fold = 2, sub = 'DESeq'), error = function(e_) { return(NA) })
  
  tryCatch(expr = dba.plotMA( dba_, contrast = 1, th = alpha_, bUsePval = F, method = DBA_EDGER,  bXY = TRUE, fold = 2, sub = 'edgeR'), error = function(e_) { return(NA) })                                        # MA plot for read distribution of DE peaks
  tryCatch(expr = dba.plotMA( dba_, contrast = 1, th = alpha_, bUsePval = F, method = DBA_DESEQ2, bXY = TRUE, fold = 2, sub = 'DESeq'), error = function(e_) { return(NA) })
  
  tryCatch(expr = dba.plotBox(dba_, contrast = 1, th = alpha_, bUsePval = F, method = DBA_EDGER,  sub = 'edgeR'), error = function(e_) { return(NA) })                       # Box plot for read distribution of DE peaks
  tryCatch(expr = dba.plotBox(dba_, contrast = 1, th = alpha_, bUsePval = F, method = DBA_DESEQ2, sub = 'DESeq'), error = function(e_) { return(NA) })
  
  graphics.off()
  
  # extracting significant DE peaks called by edgeR and DESEQ2
  message('extracting significant DE peaks called by edgeR and DESEQ2')
  
  dpea_ = rbind(as.data.frame(dba.report(dba_, method = DBA_EDGER, contrast = 1, th = alpha_, bUsePval = F)),     # edgeR is returning bullshit
                as.data.frame(dba.report(dba_, method = DBA_DESEQ2, contrast = 1, th = alpha_, bUsePval = F)))      # all significant DE peaks
  dpea_ = dpea_[!grepl(x = dpea_$seqnames, pattern = '[A-W]|Z'),]                                                   # only canonical chromosomes
  dpea_$seqnames = droplevels(dpea_$seqnames)
  dpea_$seqnames = if(!grepl(pattern = 'chr', dpea_$seqnames[1])) { paste0('chr', dpea_$seqnames) }                  # many downstream tools works when chromosomes have "chr"
  
  # merging significant DE peaks by edgeR and DESEQ2
  message('merging significant DE peaks by edgeR and DESEQ2')
  
  dup_inds = which(duplicated(dpea_[,c("seqnames","start","end")]))     # finding peak observed for a second time
  for(dup_ind in dup_inds)                                              # finds all common peaks and replaces corresponding rows in dpea_ table with the one with lowest p-value
  {
    repeat_inds = which(dpea_$seqnames %in% dpea_$seqnames[dup_ind] &     # finding both versions of repeated peaks called by edgeR and DESEQ2
                          dpea_$start %in% dpea_$start[dup_ind] &
                          dpea_$end %in% dpea_$end[dup_ind])
    low_pval_ind = repeat_inds[order(dpea_$p.value[repeat_inds])[1]]      # which one is the peak with lowest p-value?
    dpea_[repeat_inds,] = dpea_[low_pval_ind,]
  }
  dpea_ = unique(dpea_)                                                 # keeping a union (unique) of peaks across edgeR and DESEQ2
  dpea_ = dpea_[order(dpea_$Fold, dpea_$p.value),-ncol(dpea_)]          # last column is FDR which is different for edgeR and DESEQ2 and therefore inconsistent in the merged table
  
  # extracting summits of DE peaks
  
  # reading in all summit files of current contrast
  message('reading in all summit files of current contrast')
  
  fls_ = list.files(path = '../data/fastq/MACS/', pattern = '*summits.bed', full.names = T)     # all summit files of current contrast
  summits_ = list()
  for(fl_ in fls_) { summits_[[fl_]] = read.table(file = fl_, header = F, sep = '\t', col.names = c('chrom','chromStart','chromEnd','name','score')) }
  summits_ = do.call(rbind, summits_)
  summits_ = unique(summits_)                                                                                                                                                                         # keeping a union (unique) of summits across all conditions
  summits_$chrom = if(!grepl(pattern = 'chr', summits_$chrom[1])) { paste0('chr', summits_$chrom) }                  # many downstream tools works when chromosomes have "chr"
  
  # extracting summits of DE peaks
  message('extracting summits of DE peaks')
  
  summits_ls = list()       # list of summits of DE peaks
  for(row_ in 1:nrow(dpea_))
  {
    summits_ls[[row_]] = cbind(summits_[ summits_$chrom %in% dpea_$seqnames[row_] &
                                           dpea_$start[row_] <= summits_$chromStart & summits_$chromEnd <= dpea_$end[row_],],     # summit is within a peak
                               Fold = dpea_$Fold[row_])
  }
  summits_ls = do.call(rbind, summits_ls)
  summits_ls = summits_ls[order(summits_ls$Fold, summits_ls$chrom, summits_ls$chromStart),]
  
  # writing
  message('writing')
  
  write.table(dpea_, file = paste0('../out/',fl_,"_all_DEpeaks.tsv"), sep="\t", quote = F, row.names = F)                       # all significant DE peaks
  
  all_ = cbind(summits_ls[,c("chrom","chromStart","chromEnd","name","score")], strand = '+')
  write.table(all_, file = paste0('../out/',fl_,"_all_summits.bed"), sep="\t", quote = F, row.names = F, col.names = F)         # all significant DE peaks
  
  upreg_ =   cbind(summits_ls[ summits_ls$Fold < 0, c("chrom","chromStart","chromEnd","name","score")], strand = '+')             # significant up-regulated DE peaks
  write.table(upreg_,   file = paste0('../out/',fl_,"_upreg_summits.bed"), sep="\t", quote=F, row.names=F, col.names=F)         # significant up-regulated DE peaks
  
  downreg_ = cbind(summits_ls[ 0 <= summits_ls$Fold , c("chrom","chromStart","chromEnd","name","score")], strand = '+')           # significant down-regulated DE peaks
  write.table(downreg_, file = paste0('../out/',fl_,"_downreg_summits.bed"), sep="\t", quote=F, row.names = F, col.names=F)     # significant down-regulated DE peaks
}

