# Load necessary libraries
library(readr)
library(dplyr)
rm(list= ls())
# Read input files
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/4.Herbarium/4.Sequencing/8.SequencingBatch2/barcodes")
list.files()
batch2 <- read_csv("batch2_barcodes.csv")
barcode_seq <- read_csv("barcodes_sequence.csv")
lane_info <- read.csv("batch2_samples_per_lane.csv")

# Filter out samples that were prepped twice (same barcodes)
unique_samples <- batch2 %>% distinct(sample, i5_name, i7_name)

# Join i5 barcodes
batch2_i5 <- unique_samples %>%
  left_join(barcode_seq %>% filter(barcode_type == "i5") %>% select(barcode_name, i5 = sequence),
            by = c("i5_name" = "barcode_name"))

# Join i5 barcodes
batch2_i5i7 <- batch2_i5 %>%
  left_join(barcode_seq %>% filter(barcode_type == "i7") %>% select(barcode_name, i7 = sequence),
            by = c("i7_name" = "barcode_name"))

# Select and rename columns for output
output_data <- batch2_i5i7 %>%
  select(sample, i5, i7)

# Make sure no two samples have the same combination of barcodes
# Check for duplicated combinations
duplicate_combinations <- output_data %>%
  group_by(i5, i7) %>%
  filter(n() > 1)
# 16 samples with the same combination - need to keep track when demultiplexing

# Match sample barcode with lane information
output_data_lane <- merge(output_data, lane_info, by = "sample")

# Separate samples by lane - they will be demultiplexed separately
lane1 <- output_data_lane[output_data_lane$lane == "1",]
lane2 <- output_data_lane[output_data_lane$lane == "2",]


# Prep files for demultiplexing
lane1 <- read.table("lane1_barcode_sequences_with_revcomp.txt", sep = "\t", header = T)
lane2 <- read.table("lane2_barcode_sequences_with_revcomp.txt", sep = "\t", header = T)

# Export
# Write to a tab-delimited file
write_delim(lane1, "lane1_barcode_sequences.txt", delim = "\t")
write_delim(lane2, "lane2_barcode_sequences.txt", delim = "\t")


# Output format for demuxbyname.sh script
# Build the pair in header order: i7 first, then i5 !!**reverse complement**!!
lane1$pair <- paste0(lane1$i7, "+", lane1$i5_revcomp)
lane2$pair <- paste0(lane2$i7, "+", lane2$i5_revcomp)

# check
stopifnot(!any(duplicated(lane1$sample)))
stopifnot(!any(duplicated(lane1$pair)))

stopifnot(!any(duplicated(lane2$sample)))
stopifnot(!any(duplicated(lane2$pair)))


# File for demuxbyname.sh's names= argument: one barcode pair per line
writeLines(lane1$pair, "lane1_demuxbyname_batch2_barcodes.txt")
writeLines(lane2$pair, "lane2_demuxbyname_batch2_barcodes.txt")

# Lookup table to rename demuxbyname.sh's output afterward (barcode -> sample ID)
write.table(lane1[, c("pair", "sample")], "lane1_demuxbyname_batch2_barcode_to_sample.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(lane2[, c("pair", "sample")], "lane2_demuxbyname_batch2_barcode_to_sample.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

