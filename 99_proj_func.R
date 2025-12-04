dexseq_analysis <- function(isoCm, Metadata, IsoformInfo){
  n_samples <- ncol(isoCm)
  isoCm <- isoCm[rowSums(isoCm >= 10) >= (0.5 * n_samples), ]
  isoforms <- rownames(isoCm)
  IsoformInfo <- IsoformInfo |> filter(transcript_id %in% isoforms)
  IsoformInfo <- IsoformInfo |> distinct(transcript_id, .keep_all = TRUE)
  IsoformInfo <- IsoformInfo[match(rownames(isoCm),
                                   IsoformInfo$transcript_id),]
  stopifnot( all(rownames(isoCm) == IsoformInfo$transcript_id) )
  
  
  Metadata$mutation <- factor(Metadata$mutation)
  
  dxdSubset <- DEXSeqDataSet(
    countData = isoCm,
    sampleData = Metadata,
    design = ~ mutation + exon + exon:mutation,
    featureID = IsoformInfo$transcript_id,
    groupID = IsoformInfo$gene_id
  )
  
  dxdSubset <- estimateSizeFactors(dxdSubset)
  dxdSubset <- estimateDispersions(dxdSubset)
  dxdSubset <- testForDEU(
    dxdSubset,
    fullModel = ~ exon + exon:mutation,
    reducedModel = ~ exon 
  )
  dxdSubset <- estimateExonFoldChanges(dxdSubset, fitExpToVar = "mutation")
  dxdResult <- DEXSeqResults(dxdSubset)
  return(dxdResult)
}

pcoa_plot <- function(isoCm, Metadata, IsoformInfo){
  n_samples <- ncol(isoCm)
  isoCm <- isoCm[rowSums(isoCm >= 10) >= (0.5 * n_samples), ]
  isoforms <- rownames(isoCm)
  IsoformInfo <- IsoformInfo |> filter(transcript_id %in% isoforms)
  IsoformInfo <- IsoformInfo |> distinct(transcript_id, .keep_all = TRUE)
  IsoformInfo <- IsoformInfo[match(rownames(isoCm),
                                   IsoformInfo$transcript_id),]
  stopifnot( all(rownames(isoCm) == IsoformInfo$transcript_id) )
  
  
  Metadata$mutation <- factor(Metadata$mutation)
  
  dxdSubset <- DEXSeqDataSet(
    countData = isoCm,
    sampleData = Metadata,
    design = ~ mutation + exon + exon:mutation,
    featureID = IsoformInfo$transcript_id,
    groupID = IsoformInfo$gene_id
  )
  
  plotPCA(rlog(dds_ocular), intgroup=c("mutation"))
}