# Load library
library(SummarizedExperiment)
library(dplyr)

# Build SummarizedExperiment function
build_summarized_experiment <- function(combined_file,
                                        experiment_ann,
                                        lfq_type = "MaxLFQ",
                                        level = "protein"){
  # Read quant table
  prot_quant <- read.table(combined_file,
                           header = T,
                           sep = "\t",
                           quote = "",
                           comment.char = "",
                           blank.lines.skip = F,
                           check.names = F
  )
  
  # Delete contam rows
  prot_quant <- prot_quant[!grepl("contam", prot_quant$Protein),]
  
  # Take first identifier per row and make unique names
  prot_uniq <- prot_quant %>%
    mutate(
      name = if (level == "protein") get("Gene") else get("Peptide Sequence"),
      ID = if (level == "protein") get("Protein ID") else get("Peptide Sequence"),
      name = make.unique(ifelse(name == "" | is.na(name), ID, name))
      )
  
  # Set rownames
  rownames(prot_uniq) <- prot_uniq$ID
  
  # Select MaxLFQ/Intensity columns
  if (lfq_type == "MaxLFQ"){
    lfq_col <- grep("MaxLFQ", colnames(prot_uniq))
    prot_lfq <- prot_uniq[, lfq_col]
  }
  else {
    lfq_col <- grep("Intensity", colnames(prot_uniq))
    lfq_col <- lfq_col[!grepl("MaxLFQ", lfq_col)]
    prot_lfq <- prot_uniq[, lfq_col]
  }
  
  
  # Replace 0 by NA
  prot_lfq[prot_lfq == 0] <- NA
  
  # Read annotation table
  exp_anno <- read.table(experiment_ann,
                         header = T,
                         sep = "\t",
                         stringsAsFactors = F)
  if (lfq_type == "MaxLFQ"){
    exp_anno$label <- paste(exp_anno$sample, "MaxLFQ.Intensity", sep = " ")
  }
  else {
    exp_anno$label <- paste(exp_anno$sample, "Intensity", sep = " ")
  }
  
  # Set rownames
  rownames(exp_anno) <- exp_anno$label
  
  # Match column names quant with label from annotation
  matched <- match(
    make.names(exp_anno$label),
    make.names(colnames(prot_lfq))
  )
  
  # Check if labels in annotation match with column names in quant
  if (any(is.na(matched))) {
    print(make.names(exp_anno$label))
    print(make.names(colnames(prot_lfq)))
  }
  
  # Set rownames from annotation to sample name
  rownames(exp_anno) <- exp_anno$sample_name
  # Set colnames matched from quant to sample name
  colnames(prot_lfq)[matched] <- exp_anno$sample_name
  # Keep column name not NA and reorder
  prot_lfq <- prot_lfq[, !is.na(colnames(prot_lfq))][rownames(exp_anno)]
  
  # Create rowData
  row_data <- prot_uniq[, -lfq_col]
  rownames(row_data) <- prot_uniq$ID
  
  # Make SummarizedExperiment (SE)
  se <- SummarizedExperiment(
    assays = as.matrix(prot_lfq),
    colData = exp_anno,
    rowData = row_data,
    metadata = list("log2transform"=F, "lfq_type"= lfq_type,
                      "level"=level)
  )
  
  return(se)
}
