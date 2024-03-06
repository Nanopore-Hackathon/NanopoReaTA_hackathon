
args = commandArgs(trailingOnly=TRUE)
input_metadata <- args[1]
input_counttable <- args[2]
input_threads <- args[3]
output_file <- args[4]

# read metadata file (fread)
# input_metadata = path/to/metadata_file, piped in by nextflow pipeline
metadata <- fread(input_metadata)

# fread count table: input_counttable = path/to/counttable, piped in by nextflow pipeline
counts <- fread(input_counttable)

# threads as input parameter!
register(MulticoreParam(input_threads))

createDDS2 <- function(counts, metadata){
  flog.info("########## Create DDS object ###########")
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = ~ Conditions)
  
  dds <- DESeq(dds, parallel = T)
  
  return(dds)
}

run_preprocessing_dea <- function(metadata, counts, output_file){
  flog.info("########## Differential Expression Analysis ###########")
  
  # Rename metadata columns
  colnames(metadata)[c(1,2)] <- list('Samples', 'Conditions') 
  
  row.names(metadata) <- metadata$Samples
  
  missingSampleInfos = colnames(counts)[-which(colnames(counts) %in% metadata$Samples)]
  if (length(missingSampleInfos) > 0){
    print(paste0("No metadata found for the following samples: ", paste(missingSampleInfos, collapse = ",")))
    print(paste0(">>>>> Counts will be excluded!"))
    
  } else {
    print("All required information is included!")
  }
  metadata = metadata[which(row.names(metadata) %in% colnames(counts)),]
  
  metadata = metadata[which(row.names(metadata) %in% intersect(row.names(metadata), colnames(counts))),]
  
  counts = counts[, match(row.names(metadata), colnames(counts))]
  
  dds = createDDS2(counts, metadata)
  rld = rlog(dds)
  
  # save objects 
  save(dds, rld, file=output_file)
  
}

run_preprocessing_dea(input_metadata, input_counttable, output_file)

