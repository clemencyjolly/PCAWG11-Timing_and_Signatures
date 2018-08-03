# Mutational timing of gains for single samples

setwd("/home/cjolly/results/ICGC_timing") 

# libraries
library(VariantAnnotation)
library(cancerTiming)

# vcf file is input from command line 
args = commandArgs(TRUE)
vcf_filename = args[1]
id = gsub(".consensus.20160830.somatic.snv_mnv.vcf.gz","",basename(vcf_filename))
print(id)

# Assign mutations as clonal or subclonal, requires file with clustered mutations
assignMuts <- function(subclones_file){
  
  clonal.rows = NULL
  
  # define clonal cluster and subclones
  if (nrow(subclones[subclones$fraction_cancer_cells > 0.9 & subclones$fraction_cancer_cells < 1.1, ]) == 1){
    
    # if there is only one clone between 0.9 and 1.1
    # take this as the clone, and add the superclonal mutations
    # subclones are everything below 
    clonal.row = subclones[subclones$fraction_cancer_cells > 0.9 & subclones$fraction_cancer_cells < 1.1, ]
    clonal.rows = subset(subclones, fraction_cancer_cells >= clonal.row$fraction_cancer_cells)
    subclonal.rows = subset(subclones, fraction_cancer_cells < clonal.row$fraction_cancer_cells)
    no.subclones = nrow(subclonal.rows)
    
    # if there is nothing in the range of 0.9 to 1.1
    # remove superclone
    # and count only subclonal clusters and mutations 
  } else if (nrow(subclones[subclones$fraction_cancer_cells > 0.9 & subclones$fraction_cancer_cells < 1.1, ]) == 0){
    
    subclones = subset(subclones, fraction_cancer_cells < 1)
    clonal.rows = as.data.frame(matrix(ncol=3, nrow=0))
    subclonal.rows = subset(subclones, fraction_cancer_cells < 1)
    no.subclones = nrow(subclonal.rows)
    
    # if there are multiple clusters in the range of 0.9 and 1.1 
    # look between 0.95 and 1.05
  } else if (nrow(subclones[subclones$fraction_cancer_cells > 0.9 & subclones$fraction_cancer_cells < 1.1, ]) > 1){
    
    clonal.row.strict = subclones[subclones$fraction_cancer_cells > 0.95 & subclones$fraction_cancer_cells < 1.05,]
    
    # if there is only 1 cluster in the strict range 
    if (nrow(clonal.row.strict) == 1){
      
      clonal.row = clonal.row.strict
      clonal.rows = subset(subclones, fraction_cancer_cells >= clonal.row$fraction_cancer_cells)
      subclonal.rows = subset(subclones, fraction_cancer_cells < clonal.row$fraction_cancer_cells)
      no.subclones = nrow(subclonal.rows)
      
      # if there are still two clusters in the strict range
      # then take the larger ccf 
    } else if (nrow(clonal.row.strict) > 1){
      
      clonal.row = clonal.row.strict[which.max(clonal.row.strict$fraction_cancer_cells),]
      clonal.rows = subset(subclones, fraction_cancer_cells >= clonal.row$fraction_cancer_cells)
      subclonal.rows = subset(subclones, fraction_cancer_cells < clonal.row$fraction_cancer_cells)
      no.subclones = nrow(subclonal.rows)
      
      # if there is nothing between 0.95 and 1.05
      # but two within 0.9 and 1.1
    } else if (nrow(clonal.row.strict) == 0){
      
      clonal.row.relaxed = subclones[subclones$fraction_cancer_cells > 0.9 & subclones$fraction_cancer_cells < 1.1, ]
      clonal.row = clonal.row.relaxed[which.max(clonal.row.relaxed$n_snvs),]
      clonal.rows = subset(subclones, fraction_cancer_cells >= clonal.row$fraction_cancer_cells)
      subclonal.rows = subset(subclones, fraction_cancer_cells < clonal.row$fraction_cancer_cells)
      no.subclones = nrow(subclonal.rows)
    }
  }
  return(clonal.rows)
}



# dpoutput_path = path to clustering files
dpoutput_files = list.files(dpoutput_path, pattern = "cluster_assignments.txt.gz", recursive=FALSE, full.names = TRUE)
subclonal_structure_files = list.files(dpoutput_path, pattern = "subclonal_structure.txt.gz", recursive=FALSE, full.names = TRUE)

# tumour purity and ploidy
purity_ploidy = read.table("../consensus.20170217.purity.ploidy.txt.gz", header = TRUE)
sample_purity = subset(purity_ploidy, samplename == id)$purity
norm_cont = 1 - sample_purity

# cn_path = path to copy number files
cn_files = list.files(cn_path, recursive=FALSE, full.names = TRUE)  

# get mutations and read counts from vcf
sample_vcf = vcf_filename
vcf = readVcfAsVRanges(sample_vcf, "hg19",  param = ScanVcfParam(fixed=c("ALT","FILTER"),geno=NA))
print(sample_vcf)
vcf_df = as.data.frame(vcf)
vcf_df = subset(vcf_df, !is.na(vcf_df$t_alt_count) & !is.na(vcf_df$t_ref_count))
vcf_df = vcf_df[,c(1,2,3,20,21)]
  
# clustering output for specified sample
sample_dpoutput = dpoutput_files[grep(id, dpoutput_files)]  
print(sample_dpoutput)

if (file.exists(sample_dpoutput)){
  
  dp_output = read.table(gzfile(sample_dpoutput), header=TRUE)
  print("1. Sample has clustering output")
  
  # subclonal structure file for specified sample 
  sample_subclones = subclonal_structure_files[grep(id, subclonal_structure_files)]
  print(sample_subclones)
  
  if (file.exists(sample_subclones)){
    subclones = read.table(gzfile(sample_subclones), header=TRUE)
    print("2. Sample has subclonal structure file")
    
    # remove clusters with less than 1%
    subclones$n_snvs_prop = subclones$n_snvs/sum(subclones$n_snvs)  
    subclones = subset(subclones, n_snvs_prop > 0.01)   
    subclones$cluster = paste0("cluster_", subclones$cluster)
    
    clonal.rows = assignMuts(subclones)
    clonal = clonal.rows

    # remove anything not SNV and assign muts to subclones
    dp_output = subset(dp_output, mut_type == "SNV")
    clusters = unique(subclones$cluster)
    cluster.cols = subset(dp_output, select=c(clusters))
    dp_output[, "max_cluster"] = colnames(cluster.cols)[max.col(cluster.cols,ties.method="first")]
    
    clonal_mutations = subset(dp_output, max_cluster %in% clonal$cluster) 
    
    if (is.data.frame(clonal_mutations) == TRUE & nrow(clonal_mutations) > 0){
      
      print("3. Sample has clonal mutations")
      
      # subset vcf to get only clonal mutations
      vcf_gr = makeGRangesFromDataFrame(vcf_df, keep.extra.columns = T)
      clonal_mutations$start = clonal_mutations$position
      clonal_mutations$end = clonal_mutations$position
      clonal_mutations$position = NULL
      clonal_mutations_gr = makeGRangesFromDataFrame(clonal_mutations)
      
      # get overlap between vcf and clonal mutations
      hits = findOverlaps(vcf_gr, clonal_mutations_gr, type = "start")
      idx = hits@queryHits
      vcf_clonal_gr = vcf_gr[idx] 
      
      # consensus clonal copy number for specified sample
      sample_cn = cn_files[grep(id, cn_files)]
      print(sample_cn)
      
      if (file.exists(sample_cn)) {
        print("4. Sample has copy number")
        cn = read.table(gzfile(sample_cn), header=TRUE)
        
        # get clonal gain segments and reduce file 
        # take only a-d segments 
        cn_gains = subset(cn, level %in% c("a","b","c","d") & major_cn > 1)[,1:6]
        
        if (is.data.frame(cn_gains) == TRUE & nrow(cn_gains) > 0){
          print("5. Sample has clonal gains")
          
          # get mutations in regions of clonal copy number 
          cn_gains_gr = makeGRangesFromDataFrame(cn_gains, keep.extra.columns=TRUE)
          vcf_gains_gr = mergeByOverlaps(vcf_clonal_gr, cn_gains_gr, type = "within")
          muts_df = as.data.frame(vcf_gains_gr$vcf_clonal_gr)
          cn_df = as.data.frame(vcf_gains_gr$cn_gains_gr, row.names=NULL)
          
          muts_clonal_gains = cbind(muts_df, cn_df)
          names(muts_clonal_gains)[9] = "segment.start"
          names(muts_clonal_gains)[10] = "segment.end"
          
          if (is.data.frame(muts_clonal_gains) == TRUE & nrow(muts_clonal_gains) > 0){
            print("6. Sample has clonal mutations in clonal gains")
            
            # Define type of segment
            muts_clonal_gains$type = "none"
            muts_clonal_gains$type[muts_clonal_gains$major_cn == 2 & muts_clonal_gains$minor_cn == 1] = "SingleGain"
            muts_clonal_gains$type[muts_clonal_gains$major_cn == 2 & muts_clonal_gains$minor_cn == 0] = "CNLOH"
            muts_clonal_gains$type[muts_clonal_gains$major_cn == 3 & muts_clonal_gains$minor_cn == 1] = "DoubleGain"
            names(muts_clonal_gains)[1] = "chr"
            
            muts_clonal_gains$segId = paste0(muts_clonal_gains$chr,"_",muts_clonal_gains$segment.start, "_", muts_clonal_gains$segment.end, "_", muts_clonal_gains$type)
            muts_clonal_gains = subset(muts_clonal_gains, type != "none")
            
            if (is.data.frame(muts_clonal_gains) == TRUE & nrow(muts_clonal_gains) > 0){
              print("7. Sample has clonal mutations within clonal timeable gains")
              muts_clonal_gains$sample = id
              
              # do formatting for eventTiming  
              names(muts_clonal_gains)[6] = "nMutAllele"
              muts_clonal_gains$nReads = muts_clonal_gains$t_ref_count + muts_clonal_gains$nMutAllele
              muts_clonal_gains$mutationId = paste0(muts_clonal_gains$chr, "_", muts_clonal_gains$start)  
              clonal_events_list = split(muts_clonal_gains, muts_clonal_gains$sample)
              
              # timing
              arguments = list(bootstrapCI="nonparametric", minMutations = 20)
              x = eventTimingOverList(dfList = clonal_events_list, normCont = norm_cont, eventArgs = arguments)
              
              # format results
              y = getPi0Summary(x, CI=TRUE)
              piSum = na.omit(y)
              rownames(piSum) = NULL
              
              print("8. Writing output")
              
              # save both as gz files
              setwd("/home/cjolly/results/ICGC_timing/timed_segments")
              write.table(piSum, file=gzfile(paste0(id, "_timed_segments.txt.gz")), sep="\t", quote = FALSE, row.names=FALSE) 
              
              
            }
          }
        }
      }
    }
  }
} 


