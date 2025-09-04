# This script performs peak/binding site enrichment analysis using DiffBind. Its output is input to ChIPseeker

require(DiffBind)

# STEP1: Performing DEPA ####

conds_ = c('IFNy','vehicle','CD28-negative','CD28-positive')
for(cont_ in combn(length(conds_),2,simplify = F) )
{
  fl_ = paste(conds_[cont_],collapse = '_vs_')
  cat('Processing ',fl_,'\n')

  smpls_ = tryCatch(expr = suppressWarnings(read.csv(file = paste0('../data/samplesheet_',fl_,'.csv'), stringsAsFactors = T)), error = function(er) {return(NULL)})
  if(is.null(smpls_)){ next() }
  rownames(smpls_) = smpls_$SampleID

  isFIle = tryCatch(expr = suppressWarnings(load(file = paste0('../data/DiffBindObj_',fl_,'.RData'))), error = function(er) {return(NULL)})
  if(is.null(isFIle))
  {
    dba_ = dba(sampleSheet = smpls_)                                          # constructing the consensus peak set (union set) covering all peak regions in all conditions
    dba_ = dba.count(dba_, bUseSummarizeOverlaps=TRUE)                        # constructing binding affinity matrix (read count matrix)
    save(dba_, file = paste0('../data/DiffBindObj_',fl_,'.RData'))
  }

  # Extracting differentially enriched peaks ####

  peaks_dt = data.frame(dba_$binding)                                                                             # consensus peakset
  peaks_dt$CHR = paste0('chr',dba_$peaks[[1]][,"Chr"])                                                            # for unknown reason chromosome name in binding matrix is numeric!
  peaks_dt = peaks_dt[grepl(peaks_dt$CHR, pattern = '^chr([0-9]+|X|Y)', perl = T),]                                  # only canonical chromosomes
  colnames(peaks_dt)[4:5] = as.character(smpls_[colnames(peaks_dt)[4:5], "Condition"])                            # replacing column names with sample names
  peaks_dt[paste0('log2FC(',fl_,')')] = log2(peaks_dt[4]/peaks_dt[5])                                             # log2FC
  peaks_dt = peaks_dt[1 <= abs(peaks_dt[,6]),]                                                                    # 2 fold difference is considered significantly differentially enriched
  peaks_dt$summit = peaks_dt$START + dba_$summits                                                                 # summit
  peaks_dt$name = apply(X = peaks_dt[4:5], MARGIN = 1, FUN = function(rw_){ return(names(which.max(rw_))) })      # name
  peaks_dt$strand = '+'                                                                                           # strand
  
  # Extracting up-/down-regulated peaks (with increased affinity/read count) ####
  
  upreg_dt = peaks_dt[0 < peaks_dt[,6],]
  downreg_dt = peaks_dt[peaks_dt[,6] <= 0,]
  
  # visualization ####
  labels_ = c(Upregulated   = round(nrow(upreg_dt)/nrow(peaks_dt),4)*100,
              Downregulated = round(nrow(downreg_dt)/nrow(peaks_dt),4)*100)
  
  pdf(file = paste0('../out/',fl_,".pdf"), bg = 'white')
  pie(labels_, labels = paste0(labels_,'%'), main = paste(conds_[cont_],collapse = ' vs '), col = rainbow(length(labels_)))
  legend("topright", names(labels_), cex = 0.8,fill = rainbow(length(labels_)))
  graphics.off()
  
  # Writing ####
  
  write.table(peaks_dt[,c(1:3,8,6,9,4,5,7)], file = paste0('../out/',fl_,".tsv"), sep="\t", quote = F, row.names = F, col.names = T)     # all significant DE peaks
  write.table(upreg_dt[,c(1:3,8,6,9,4,5,7)], file = paste0('../out/',fl_,"_upreg.tsv"), sep="\t", quote = F, row.names = F, col.names = T)     # all significant DE peaks
  write.table(downreg_dt[,c(1:3,8,6,9,4,5,7)], file = paste0('../out/',fl_,"_downreg.tsv"), sep="\t", quote = F, row.names = F, col.names = T)     # all significant DE peaks
}

