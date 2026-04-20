install.packages(c("shiny", "ggplot2", "dplyr", "tidyr", "readr"))
BiocManager::install("org.Hs.eg.db", update = FALSE)


suppressPackageStartupMessages({
  library(shiny) 
  library(ggplot2) 
  library(dplyr) 
  library(tidyr) 
  library(readr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# 1. Data Reading and Cleaning
data <- read_tsv("GSE104296_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")

data_df <- as.data.frame(data)
rownames(data_df) <- data_df$GeneID
data_clean <- data_df[, -1]

# 2. Create Metadata
meta <- data.frame(
  sample = c("GSM2794663", "GSM2794664", "GSM2794665", "GSM2794666", 
             "GSM2794667", "GSM2794668", "GSM2794669", "GSM2794670"),
  condition = c("Control", "Control", "ERb-KO", "ERb-KO", 
                "ERb1", "ERb1", "ERb5", "ERb5")
)

# 3. Gene Annotation
all_ids <- as.character(rownames(data_clean))
mapping_results <- AnnotationDbi::select(org.Hs.eg.db, 
                                         keys = all_ids, 
                                         columns = c("SYMBOL", "GENENAME"), 
                                         keytype = "ENTREZID")

# 4. Processing results
mapping_table <- mapping_results %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
  dplyr::rename(GeneID = ENTREZID, 
                Symbol = SYMBOL, 
                `Gene Name / Type` = GENENAME)

# 5. Calculating Statistics
control_cols <- c("GSM2794663", "GSM2794664")
ko_cols      <- c("GSM2794667", "GSM2794668")

# Calculating Log2 Fold Change and p-values
stats_df <- data_clean %>%
  mutate(
    log2FC = rowMeans(across(all_of(ko_cols))) - rowMeans(across(all_of(control_cols))),
    # Performing t-test
    pvalue = apply(as.matrix(.), 1, function(x) {
      # Handling genes with zero variance/constant values
      tryCatch(t.test(x[ko_cols], x[control_cols])$p.value, error = function(e) NA)
    })
  ) %>%
  tibble::rownames_to_column("GeneID") %>%
  dplyr::select(GeneID, log2FC, pvalue)

final_mapping_table <- mapping_table %>%
  left_join(stats_df, by = "GeneID")

write_csv(final_mapping_table, "annotated_with_stats.csv")

# Saving cleaned expression data for the app
write_csv(tibble::rownames_to_column(data_clean, "GeneID"), "expression_data_clean.csv")