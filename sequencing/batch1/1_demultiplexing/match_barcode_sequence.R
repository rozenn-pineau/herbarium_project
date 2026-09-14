# Load necessary libraries
library(readr)
library(dplyr)
rm(list= ls())
# Read input files
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/4.Herbarium/4.Sequencing/6.SequencingBatch1/barcodes")
batch1 <- read_csv("batch1_barcodes.csv")
barcode_seq <- read_csv("barcode_sequence_correspondence.csv")

# Filter out samples that were prepped twice (same barcodes)
unique_samples <- batch1 %>% distinct(sample, i5_name, i7_name)

# Join i5 barcodes
batch1_i5 <- unique_samples %>%
  left_join(barcode_seq %>% filter(barcode_type == "i5") %>% select(barcode_name, i5 = sequence),
    by = c("i5_name" = "barcode_name"))

# Join i5 barcodes
batch1_i5i7 <- batch1_i5 %>%
  left_join(barcode_seq %>% filter(barcode_type == "i7") %>% select(barcode_name, i7 = sequence),
            by = c("i7_name" = "barcode_name"))

# Select and rename columns for output
output_data <- batch1_i5i7 %>%
  select(sample, i5, i7)

# Make sure no two samples have the same combination of barcodes
# Check for duplicated combinations
duplicate_combinations <- output_data %>%
  group_by(i5, i7) %>%
  filter(n() > 1)
# empty file - ouf


# output formats
# Write to a tab-delimited file
write_delim(output_data, "output_barcode_sequences.txt", delim = "\t")

# Create i5 FASTA file
i5_fasta <- paste0(">", output_data$sample, "\n", output_data$i5)
writeLines(i5_fasta, "i5_barcodes.fasta")

# Create i7 FASTA file
i7_fasta <- paste0(">", output_data$sample, "\n", output_data$i7)
writeLines(i7_fasta, "i7_barcodes.fasta")

# Output format for demuxbyname.sh script
# Build the pair in header order: i7 first, then i5
output_data$pair <- paste0(output_data$i7, "+", output_data$i5)

stopifnot(!any(duplicated(output_data$sample)))
stopifnot(!any(duplicated(output_data$pair)))

# File for demuxbyname.sh's names= argument: one barcode pair per line
writeLines(output_data$pair, "demuxbyname_batch1_barcodes.txt")

# Lookup table to rename demuxbyname.sh's output afterward (barcode -> sample ID)
write.table(output_data[, c("pair", "sample")], "demuxbyname_batch1_barcode_to_sample.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

