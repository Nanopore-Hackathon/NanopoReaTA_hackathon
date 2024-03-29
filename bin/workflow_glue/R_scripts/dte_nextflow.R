#_______________________________________________________________________________
# CREATE DRIMSEQ OBJECT IN NEXTFLOW

args <- commandArgs(trailingOnly = TRUE)

metadata_path = args[1] # First column Sample names, second column condition
counts_path = args[2] # counts file from transcripts 
gtf_file = args[3] # converted_gtf.csv
output_RDATA = args[4] # output path to DTU_objects.RData object

process_input_files <- function(counts_path, metadata_path){
  counts_df = data.table::fread(counts_path, data.table = T)
  geneNames = counts_df$Name
  counts_df = counts_df %>% dplyr::select(-1:-4)
  counts_df = as.data.frame(counts_df)
  counts_df_int = sapply(counts_df, as.numeric)
  row.names(counts_df_int) = geneNames
  counts_df_int <- as.data.frame(counts_df_int)
  
  samps = data.table::fread(metadata_path, data.table = F)
  samps["Samples"] = as.character(samps[,1])
  samps["Condition"] = as.character(samps[,2])
  return(list(samps, counts_df_int))
}

gtf_file = data.table::fread(gtf_path, data.table = T)


create_DRIMSeq_object <- function(counts_df, samps, gtf_file, output_RDATA){
  
  ###############################################################################################
  #                                                                                             #
  #                                                                                             #
  #                                                                                             #
  #                             Read count file, gtf and metadata                               #
  #                                                                                             #
  #                                                                                             #
  #                                                                                             #
  #                                                                                             #
  ###############################################################################################

  counts_df = counts_df[,which(colnames(counts_df) %in% samps$Samples)]

  all(rownames(counts_df) %in% gtf_file$transcript_id)
  gtf_file <- gtf_file[match(rownames(counts_df),gtf_file$transcript_id),]
  all(rownames(counts_df) == gtf_file$transcript_id)
  
  counts <- data.frame(gene_id=gtf_file$gene_id,
                       feature_id=gtf_file$transcript_id,
                       counts_df)    
  d <- dmDSdata(counts=counts, samples=samps)
  n <- length(samps$Samples)
  n.small <- min(table(samps$Condition))
  out2 <- tryCatch(
    {d <- dmFilter(d,
                   min_samps_feature_expr=as.integer(n.small), min_feature_expr=5,
                   min_samps_gene_expr=(n.small), min_gene_expr=20)
    x = F
    },
    error = function(e){
      x = T
    },
    finally = {
    })
  if (out2){
    flog.info("########## DrimSeq Filtering failed ###########")
    flog.info("Either only one splicing variant for every gene in dataset or it must be sequenced deeper")
    return(NULL)
  }
  
  input_design <- DRIMSeq::samples(d)
  return(d)
}


createDEXSeq_object = function(d_list){
  output_list = list()
  d = d_list$drim
  
  sample.data <- DRIMSeq::samples(d)
  print(sample.data)
  sample.data$Condition <- factor(sample.data$Condition,levels = c(first.level,ref.level))
  count.data <- round(as.matrix(counts(d)[,-c(1:2)]))

  dxd <- DEXSeqDataSet(countData=count.data,
                       sampleData=sample.data,
                       design=~Samples + exon + Condition:exon,
                       featureID=counts(d)$feature_id,
                       groupID=counts(d)$gene_id)
  dxd$Condtition = factor(dxd$Condition,levels = c(first.level,ref.level))
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet=TRUE)
  dxd <- testForDEU(dxd, reducedModel=~Samples + exon)
  return(dxd)
}


meta_counts = process_input_files(counts_path, metadata_path)

d_object = create_DRIMSeq_object(counts_df = meta_counts[2], samps = meta_counts[1], gtf_file, output_RDATA)
dxd_object = createDEXSeq_object(d_list = d_object)

save(d_object, dxd_object, file = output_RDATA)

