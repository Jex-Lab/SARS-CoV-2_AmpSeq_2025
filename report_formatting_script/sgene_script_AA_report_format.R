
library(tibble)
library(writexl)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(lubridate)
library(dplyr)

# -----------------------------------------------------------------
# Summary table
# -----------------------------------------------------------------

# import Sample List file exported from Benchling
sample_info <- read.csv("Sample_List.csv")[, c(2,7,8,10,16,17)]

rownames(sample_info) <- NULL

# set column names
colnames(sample_info) <- c("RNA_ID",
                           "Sample ID",
                           "DHHS SiteName",
                           "Retrieved",
                           "Sample Type",
                           "Sample Units")

# -----------------------------------------------------------------

# import read_counts
read_counts <- read.table("read_counts.txt", sep = "\t", header = TRUE)

# ensure filtered read values are >0 so percentage calculations don't fail
read_counts$filtered <- read_counts$filtered + 1

# calculate percentages
read_counts[,c(3:ncol(read_counts))] <- round((read_counts[,c(3:ncol(read_counts))]/read_counts$filtered)*100, digits = 1)

# calculate Unassigned
read_counts$Unassigned <- 100 - rowSums(read_counts[,c(3:(ncol(read_counts)))])

# set column names
colnames(read_counts)[1] <- "RNA_ID"
colnames(read_counts)[2] <- "Filtered Reads"

# -----------------------------------------------------------------

# import snp_counts - column names will be offset
snp_counts <- read.table("snp_counts.txt", sep = "\t", header = TRUE, fill = TRUE, row.names = NULL)
snp_counts <- snp_counts[,1:4]

# set column names
colnames(snp_counts) <- c("RNA_ID",
                          "Non-classified SNPs",
                          "Non-classified SNPs %",
                          "Detected SNPs")

# -----------------------------------------------------------------

# combine data
all_amp <- merge(sample_info, (merge(read_counts, snp_counts,by="RNA_ID")),by = "RNA_ID")

# remove NTC samples
all_amp <- all_amp[grep("NTC", all_amp$RNA_ID, invert = TRUE),]

# remove Non-Reportable samples
all_amp <- all_amp[grep("nr", all_amp$RNA_ID, invert = TRUE),]

# reset row numbers
rownames(all_amp) <- NULL

# -----------------------------------------------------------------

# set default technical column values
`Technical Notes` <- NA
`Sequencing QC` <- "Fail"

# add columns
all_amp <- cbind(all_amp, `Technical Notes`, `Sequencing QC`)

# insert technical note
all_amp$`Technical notes`[all_amp$`Filtered Reads` < 2000] <- "QC criteria not met (Filtered Reads ≥ 2k)"

# insert sequencing QC pass/fail comment
all_amp$`Sequencing QC`[is.na(all_amp$`Technical notes`)] <- "Pass"

# insert default Repeat Number
all_amp <- add_column(all_amp, `Repeat Number` = NA, .after = 1)

# -----------------------------------------------------------------

# reorder columns
all_amp <- all_amp[,c("RNA_ID",
                      "Repeat Number",
                      "Sample ID",
                      "DHHS SiteName",
                      "Retrieved",
                      "Sample Type",
                      "Sample Units",
                      "Filtered Reads",
                      "Sequencing QC",
                      "NB181",
                      "XEC",
                      "LF7",
                      "MV1",
                      "LP81",
                      "XFG",
                      "PY1",
                      "PA1",
                      "NY6",
                      "PQ6",
                      "NY10",
                      "Unassigned",
                      "Non-classified SNPs",
                      "Non-classified SNPs %",
                      "Detected SNPs",
                      "Technical notes")]

# -----------------------------------------------------------------
# Identify the 10 most abundant "Unassigned"
# -----------------------------------------------------------------

# import other count data
other_counts <- read.table(file = "other_counts.txt", header = TRUE)
colnames(other_counts)[1] <- c("RNA_ID")
other_counts_only <- ncol(other_counts)

# retain only the reportable samples
other_counts <- merge(other_counts, all_amp, by = "RNA_ID")
other_counts <- other_counts[,1:other_counts_only]

# get identities of reads described in top10_other_snps.txt from
# https://cov-spectrum.org/explore/World/AllSamples/AllTimes/
colnames(other_counts)[3] <- c("XEC -T478K")
colnames(other_counts)[4] <- c("NB.1.8.1 L452W>R")
colnames(other_counts)[5] <- c("LP.8.1.1 -T478K")
colnames(other_counts)[6] <- c("LP.8.1.1 T478K>I")
colnames(other_counts)[7] <- c("XEC +N487D")
colnames(other_counts)[8] <- c("XFG T478K>I")
colnames(other_counts)[9] <- c("XFG -K444R")
colnames(other_counts)[10] <- c("PY.1 -L441R")
colnames(other_counts)[11] <- c("XFG -N487D")
colnames(other_counts)[12] <- c("XFG F456L>R")

# calculate percentages
other_counts[,c(3:ncol(other_counts))] <- round((other_counts[,c(3:ncol(other_counts))]/other_counts$filtered)*100, digits = 1)

# calculate percentages of 'Other' - the percentage not accounted for by the binned count data or the identified Unassigned count data
read_counts_only <- ncol(read_counts)
read_counts <- merge(read_counts, all_amp, by = "RNA_ID")
read_counts <- read_counts[,1:read_counts_only]
other_counts$Other <- read_counts$Unassigned - rowSums(other_counts[,3:ncol(other_counts)])

# -----------------------------------------------------------------

# shorten sample names in the final data-frames
all_amp$RNA_ID <- substr(all_amp$RNA_ID,6,11)
other_counts$RNA_ID <- substr(other_counts$RNA_ID,6,11)

# export data
write_xlsx((merge((all_amp[,c(1,2,9,21)]), (other_counts[,c(1,3:ncol(other_counts))]), by = "RNA_ID")),
           "unassigned_sheet.xlsx")

write_xlsx(all_amp,
           "summary_sheet.xlsx")

# -----------------------------------------------------------------
# "Unassigned" bar chart
# -----------------------------------------------------------------

# retain only PASS samples
other_counts <- other_counts[other_counts$filtered >= 2000,]

# convert data to long format (exclude "filtered" column)
other_long <- reshape2::melt(other_counts[,-2])

# # legend title
# legend_title <- "Description"

# create colour palette
nb.cols <- length(c(3:ncol(other_counts))) # number of colours needed
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

other_long %>%
  ggplot(aes(fill=variable, y=value, x=RNA_ID)) +
  geom_bar(position="stack",stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  ylab("Relative Abundance (%)") +
  xlab("Sample") +
  scale_fill_manual(name=NULL, values=mycolors)
  # scale_fill_manual(legend_title, values=mycolors)

ggsave("Unassigned_BarChart.pdf", width = 420, height = 210, units = "mm")

# -----------------------------------------------------------------
