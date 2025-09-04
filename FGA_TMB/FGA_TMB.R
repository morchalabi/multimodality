# This script computes FGA and TMB using ichorCNA's tumor.seg.txt files and lpWGS bam files, respectively.
# TMB was computed on variants called by GATK Mutect2: https://gatk.broadinstitute.org/hc/en-us/articles/30331989211419--Tool-Documentation-Index
# For how the annotation of variants, output by FilterMutectCalls, is done by Funcotator read 'Output' section from this
# page: https://gatk.broadinstitute.org/hc/en-us/articles/30332018805787-Funcotator

library(httr)     # for Genome Table POST

# Fraction of genome altered (FGA) ####

fga_ = numeric()
patients_ = conditions_ = NULL
for(f_ in list.files(path = '../data', pattern = '(tumor.seg.txt)', recursive = T, full.names = T))
{
  # patient IDs
  ptn_ = regexpr(pattern = 'W[0-9]+_', text = f_, perl = T)
  ptn_ = paste(strsplit(f_, split = '')[[1]][ptn_[1]: (ptn_[1] + attr(ptn_,"match.length")-2)], collapse = '')
  patients_ = c(patients_, ptn_)
  
  # treatment condition
  cond_ = regexpr(pattern = '(neoadjuvant)|(untreated)', text = f_, perl = T)
  cond_ = paste(strsplit(f_, split = '')[[1]][cond_[1]: (cond_[1] + attr(cond_,"match.length")-1)], collapse = '')
  conditions_ = c(conditions_, cond_)
  
  # load ichorCNA segments file
  cna_segments = read.delim(file = f_, header = TRUE, sep = "\t", check.names = T, as.is = T, quote = "")

  # filter CNA regions to exclude neutral copy number (Corrected_Copy_Number == 2)
  filtered_cnas = subset(cna_segments, Corrected_Call != 'NEUT')

  # calculate the total length of altered regions
  cna_length = sum(filtered_cnas$end - filtered_cnas$start)

  # total mappable genome length (~3 billion bases)
  total_genome_length = 3e9

  # compute FGA
  fga_ = c(fga_, round(cna_length / total_genome_length, digits = 2)*100)
}

# Tumor mutational burden (TMB) ####

tmb_ = numeric()                          # tumor mutational burden
somatic_mutations_ls = list(); j_ = 1     # a list storing genes affected by somatic mutations of all patients
for(f_ in list.files(path = '../data', pattern = '(filtered_PASS-weak-evidence_annotated.vcf.gz)', recursive = T, full.names = T))
{
  cat('reading vcf file of ',f_, '\n')
  
  # extracting non-synonymous mutations (needed for TMB) of current patient
  
  lines_ = readLines(con = f_)            # reading lines of vcf files output by GATK's FilterMutectCalls
  somatic_mutations = list(); i_ = 1      # a list storing genes affected by somatic mutations of current patient
  for(l_ in lines_)                       # for every line in the vcf file
  {
    # extracting annotation related to gene mutation ####
    
    annot_ = regexpr(text = l_, pattern = '(FUNCOTATION=\\[).+(g.chr)', perl = T)     # part of the annotation line added by GATK's Funcotator
    if(annot_ == -1) next()     # no match?
    
    length_ = attr(annot_,"match.length")-1
    
    l_1 = strsplit(x = l_, split = '')[[1]]
    l_1 = paste(l_1[annot_: (annot_+length_)], collapse = '')
    l_1 = strsplit(x = l_1, split = '\\[|\\|')[[1]][c(2,4:6,7:9,11:12)]
    gene_ = l_1[1]                # gene affected
    chr_ = l_1[2]                 # chromosome
    str_ = as.integer(l_1[3])     # start
    end_ = as.integer(l_1[4])     # end
    var_class = l_1[5:6]          # class of variant; if var_class[1] is SPLICE_SITE, there could be additional variant classification like 'SPLICE_SITE,INTRON'
    var_type = l_1[7]             # type of variant
    ref_ale = l_1[8]              # reference allele for the position at which this this variant allele occurs. For insertions, this will be set to "-".
    alt_ale = l_1[9]              # alternative allele; For deletions, this will be set to "-".
    
    # extracting other stuff ####
    
    if(any(var_class %in% c('MISSENSE','NONSENSE','NONSTOP','IN_FRAME_DEL','IN_FRAME_INS','FRAME_SHIFT_INS','FRAME_SHIFT_DEL','START_CODON_SNP','START_CODON_INS','START_CODON_DEL','SPLICE_SITE')) &
       !var_class[2] %in% c('SILENT','INTRON'))
    {
      # extracting annotation related to TLOD, AD, AF and DP
      
      annot_ = regexpr(text = l_, pattern = '(TLOD=).+', perl = T)      # part of the annotation line added by GATK's Funcotator
      length_ = attr(annot_,"match.length")-1
      
      l_2 = strsplit(x = l_, split = '')[[1]]
      l_2 = paste(l_2[annot_: (annot_+length_)], collapse = '')
      l_2 = strsplit(x = l_2, split = '\t')[[1]]
      
      tlod_ = as.numeric(strsplit(x = l_2[1], split = '=')[[1]][2])     # Tumor Log of Odds (LOD): log-10 likelihood ratio that the variant is a tumor variant and not a germline one
      metrics_ = strsplit(x = l_2[2], split = ':')[[1]]                 # metrics of AD, AF and DP
      ids_ = which(metrics_ %in% c("AD","AF","DP"))
      metrics_ = strsplit(x = l_2[3], split = ':')[[1]][ids_]           # which indices?
      ad_ = metrics_[1]                                                 # Allele Depth: number of reads supporting the reference allele and the alternate allele, respectively
      af_ = as.numeric(metrics_[2])                                     # Allele Frequency: fraction of reads supporting the alternate allele over the total reads.
      dp_ = as.integer(metrics_[3])                                     # Depth of Coverage: total number of reads covering this position.
      
      # This filtering is necessary for lpWGS data ####
      # if(dp_ < 3 | tlod_ < 8 | af_ < 0.7) next()
      if(dp_ < 3 | tlod_ < 7 | af_ < 0.7) next()
      
      ############### REMOVE #############
      # if(gene_ == 'FMO2')
      # {
      #   cat('hi')
      # }else
      # {
      #   next()
      # }
      ####################################
      
      # Broad's reference vcf file (af-only-gnomad.hg38.vcf.gz) has removed all nonsynonymous alleles (variants) as pathogenic, while original gnomAD still retains many of
      # them as well-tolerated ones! This led to an excessive number of tumor alleles called by Mutect2/FilterMutectCalls. Below we remove tumor alleles well tolerated by population ####
      
      ref_AF = 0                            # if current tumor allele exists as an alternative allele at reference gnomAD, what is its population frequency (to gauge its tolerability)? Otherwise 0
      if(!var_type %in% c('INS','DEL'))     # an alternative allele cannot be of type indels and yet well-tolerated by population. gnomAD contains a few indels but with very low frequency!
      {
        # send the POST request to UCSC table browser to extract reference AF of current alternative allele (variant)
        response_ = NULL
        query_ = 'FAILED'
        while(query_ != "Success")                                                                            # Check for success; the server may be busy
        {
          message('\n\t<< contacting USCS Table Browser ... >>')
  
          response_ = POST(# UCSC Table Browser URL
                          url = "http://genome.ucsc.edu/cgi-bin/hgTables",
                          # query parameters
                          body = list(db = "hg38",                                    # genome build
                                      hgta_group = "varRep",                          # group: Variation and Repeats
                                      hgta_track = "gnomadGenomesVariantsV4_1",       # track: gnomAD Genomes Variants v4.1
                                      hgta_table = "gnomadGenomesVariantsV4_1",       # table
                                      hgta_regionType = "range",                      # testrict to a specific range
                                      position = paste0(chr_,':',str_,'-',end_),      # genomic range
                                      hgta_outputType = "primaryTable",               # output as a primary table
                                      boolshad.hgta_doTopSubmit = "get output"),      # submit the query
                          encode = "form")
          query_ = http_status(response_)$category
          if(query_ != "Success"){ cat('\ttrying again'); Sys.sleep(0.25) }
        }
        result_text = content(response_, "text")                                                              # extract the content as text
        result_table = read.delim(text = result_text, sep = "\t", header = TRUE, quote = "", as.is = T,
                                  check.names = F, colClasses = c(ref = 'character', alt = 'character'))      # convert to a data frame
        
        result_table = result_table[(str_ <= result_table$thickEnd & result_table$thickEnd <= end_),]         # each base in UCSC GT/GB is represented by an interval (thick) with [str-1,str]. Also, POST returns every thicks
                                                                                                              # containing str. This line makes sure only thicks falling within this current locus query ()
                                                                                                              # even if result_table is already empty, this line returns empty table with no error

        if(0 < nrow(result_table))                                                                            # result_table may have several alleles unrelated to current tumor one given current genomic range by chr_, st_ and end_
        {
          # tumor allele is an SNP like chr13:38848492-38848492. However, this locus on reference genome may include several alternative alleles.
          # For example result_table may have 3 rows of ref an alt alleles:
          # T G
          # T A -> this what we look for!
          # T c
          if(end_ == str_)
          {
            result_table = result_table[ result_table$alt %in% alt_ale,]      # only the record having current tumor allele
            if(0 < nrow(result_table)){ ref_AF = result_table$AF }            # frequency of tumor alternative allele in the reference population
          }
          
          # if tumor allele is a multi-nucleotide polymorphism (DNP, TNP or MNP; see variantType of GATK Funcotator), there is one AF for each base (SNP), e.g. chr2:231222761-231222762 or chr17:45241411-45241412 on GB.
          # For example result_table may return 5 rows of (ref, alt) pairs given current locus for a DNP like (ref:AC, alt:TG):
          # A C
          # C G -> this is what we need
          # A T -> this is what we need
          # C T
          # A A
          if(str_ < end_)
          {
            # creating pairs of (ref_i,alt_i) as columns of a table just like POST(). For example, for (ref:TCG, alt:ACC) it returns: (T,A), (C,C) and (G,C)
            ref_ale_alt_alt = data.frame(ref = strsplit(ref_ale, split = '')[[1]],
                                         alt = strsplit(alt_ale, split = '')[[1]])
            ref_ale_alt_alt$str = str_:(str_+nrow(ref_ale_alt_alt)-1)     # locus of each pair; for an MNP like TT -> CC result_table will return (T,C) and (T,C); these are separable only by locus (chr1:23875429-23875430)
            
            ref_AF_tmp = NULL                                             # frequency of alt_i in the reference population
            for(r_ in 1:nrow(ref_ale_alt_alt))                            # for each pair
            {
              # from result table, it returns the row matching (ref_i,alt_i)
              tmp_ = result_table[result_table$ref %in% ref_ale_alt_alt$ref[r_] &           # ref_i
                                  result_table$alt %in% ref_ale_alt_alt$alt[r_] &           # alt_i
                                  result_table$thickEnd %in% ref_ale_alt_alt$str[r_],]      # for an MNP like TT -> CC result_table will return (T,C) and (T,C); these are separable only by locus (chr1:23875429-23875430)
              ref_AF_tmp[r_] = 0
              if(nrow(tmp_) == 1)                                                           # if alt_i exists in gnomAD
              {
                # gnomAD has variant classification per base, while Funcotator assign each MNP only one classification. For example, if a DNP is (missense, synonymous) on gnomAD, Funcotator returns only missense.
                # Check this out: chr1:171199445-171199446 for W612.
                if(!tmp_$annot %in% 'synonymous'){ ref_AF_tmp[r_] = tmp_$AF }else{ ref_AF_tmp[r_] = 1 }     # if alt_i is synonymous then its AF is set to 1 to be ignored
              }
            }# for each (ref_i,alt_i)
            ref_AF = prod(ref_AF_tmp)                                     # the probability of observing an alternative MNP is the product of each base's frequency (alt_i)
          }
        }
      }

      # putting everything together
      somatic_mutations[[i_]] = data.frame(gene = gene_, chr = chr_, str = str_, end = end_,
                                           ref_allele = ref_ale, alt_allele = alt_ale,
                                           variant_class = paste(var_class[var_class != ''], collapse = ','),
                                           variant_type = var_type,
                                           TLOD = tlod_, AD = ad_, AF = af_, DP = dp_,
                                           ref_AF = ref_AF)
      i_ = i_ + 1
      cat(i_,'\n\n')
    }# nonsynonymous variant
  }# currebt patient
  somatic_mutations = do.call(somatic_mutations, what = rbind)
  
  # calculating TMB
  tmb_ = c(tmb_, round(nrow(somatic_mutations[somatic_mutations$ref_AF <= (1/1000),])/30, 2))     # an alternative variant with a AF greater than 1/1000 (0.001), by convention, is not deemed to be a tumor allele but a reference alternative variant!
                                                                                                  # total length of coding regions (exons) in the human genome is approximately 30 million base pairs (Mbp).
                                                                                                  # The UCSC Genome Browser's gnomadGenomesVariantsV4_1 comprises variant data from 76,215 WGS
  # adding current patient to the list of all patients
  somatic_mutations_ls[[f_]] = data.frame(patient = patients_[j_], condition = conditions_[j_], somatic_mutations)
  j_ = j_ + 1
}

# Writing ####

somatic_mutations_ls = do.call(somatic_mutations_ls, what = rbind)
write.table(x = somatic_mutations_ls, file = '../out/somatic_mutation_genes.tsv',sep = '\t', quote = F, row.names = F, col.names = T)

dt_ = data.frame(Patient = patients_, condition = conditions_, FGA = fga_, TMB = tmb_)
write.table(x = dt_, file = '../out/FGA_TMB.tsv',sep = '\t', quote = F, row.names = F, col.names = T)

