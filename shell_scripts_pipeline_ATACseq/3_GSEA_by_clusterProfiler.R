# This script performs gene set enrichment analysis regulated by selected peaks using clusterProfiler.

library(clusterProfiler)
require(org.Hs.eg.db)
library(enrichplot)

i_ = 1; p = list()
fls_ = list.files(path = '../out', pattern = '_annot', full.names = F)
for(fl_ in fls_)
{
  message(fl_)
  
  # STEP 1: Reading in annotation file ####
  
  annot_ = read.table(file = paste0('../out/', fl_), header = T, sep = '\t', as.is = T, check.names = T, quote = '')
  
  # STEP 2: Run GO enrichment analysis ####
  
  ego = enrichGO(gene = annot_$geneId,
                 keyType = "ENTREZID",
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 pvalueCutoff = .2,
                 readable = TRUE)
  
  ekegg = enrichKEGG(gene = annot_$geneId,
                     organism = 'hsa',
                     pvalueCutoff = .2)

  # STEP 3: Dotplot visualization ####
  
  pdf(file = paste0("../out/",strsplit(fl_, split = '[.]')[[1]][1],"_GSEA.pdf"), width = 15, height = 10)
  p = dotplot(ego, showCategory = 50)
  tryCatch(expr = suppressWarnings(plot(p, main = 'GO')), error = function(err){ return(NULL)})
  p = dotplot(ekegg, showCategory = 50)
  tryCatch(expr = suppressWarnings(plot(p, main = 'KEGG')), error = function(err){ return(NULL)})
  graphics.off()
  
  # STEP 4: Output results from GO analysis to a table ####
  
  cluster_summary = data.frame(ego)
  write.table(cluster_summary, file = paste0("../out/",strsplit(fl_, split = '[.]')[[1]][1],"_GSEA.tsv"), sep = '\t', quote = F, row.names = F, col.names = T)
}


