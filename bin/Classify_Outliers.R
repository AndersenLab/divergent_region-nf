#!/usr/bin/env Rscript
library(tidyverse)
library(anomalize)
library(tibbletime)

# args 
# 1 - sample name
# 2 - genomic region bed file
# test: args <- c("AB1","ARMS_CENTERS.bed.gz")
args <- commandArgs(trailingOnly = TRUE)

# load start and stop information for chromosomal regions
genomic_regions <- readr::read_tsv(args[2], col_names = c("CHROM", "START_REGION", "END_REGION", "REGION_TYPE"))

# load variant counts
count_df <- data.table::fread(glue::glue("{args[1]}_variant_counts.txt"), col.names = c("CHROM", "START_BIN", "END_BIN", "COUNT")) %>%
  dplyr::mutate(MID_BIN = ((END_BIN - START_BIN)/2)+START_BIN) 

# define centers, arms, and telomeric regions
for(chrom in unique(genomic_regions$CHROM)){
  # define regions
  temp_Larm <- dplyr::filter(genomic_regions, CHROM == chrom) %>%
    dplyr::filter(REGION_TYPE == "ARM", START_REGION == min(START_REGION))
  temp_Rarm <- dplyr::filter(genomic_regions, CHROM == chrom) %>%
    dplyr::filter(REGION_TYPE == "ARM", START_REGION == max(START_REGION))
  temp_center <- dplyr::filter(genomic_regions, CHROM == chrom, REGION_TYPE == "CENTER")
  
  # filter wholegenome dataset to chrom
  chrom_counts <- count_df %>%
    dplyr::filter(CHROM == chrom) %>%
    dplyr::ungroup()
  
  # make fake date to use tibbletime and anomalize functionality
  fake_time <- tibbletime::create_series('1900' ~ '2017', 'daily', class = "Date")[1:nrow(chrom_counts),]
  
  # find outliers
  temp_outlier_time <- chrom_counts %>%
    dplyr::mutate(date = fake_time$date) %>%
    as.tbl() %>%
    time_decompose(COUNT, method = "stl", trend = 100) %>%
    anomalize(remainder, method = "gesd", alpha = 0.005, max_anoms = 0.05, verbose = T)
  
  # extract outlier information
  outlier_df <- temp_outlier_time$anomaly_details$outlier_report %>%
    dplyr::select(index, rank, outlier, direction)
  
  # append outlier information count df
  outlier_by_chrom <- chrom_counts %>%
    dplyr::ungroup() %>%
    dplyr::mutate(index = 1:n()) %>%
    dplyr::left_join(.,outlier_df, by = "index")
  
  # deal with NAs from join
  outlier_by_chrom$outlier[is.na(outlier_by_chrom$outlier)] <- "No"
  outlier_by_chrom$direction[is.na(outlier_by_chrom$direction)] <- "NA"
  
  # append chromosomal region information
  if(!exists("count_regions")){
    count_regions <- outlier_by_chrom %>%
      dplyr::rowwise() %>%
      dplyr::mutate(GENOMIC_REGION = ifelse(CHROM == chrom && MID_BIN > temp_Larm$START_REGION && MID_BIN < temp_Larm$END_REGION, "LEFT_ARM",
                                            ifelse(CHROM == chrom && MID_BIN > temp_Rarm$START_REGION && MID_BIN < temp_Rarm$END_REGION, "RIGHT_ARM",
                                                   ifelse(CHROM == chrom && MID_BIN > temp_center$START_REGION && MID_BIN < temp_center$END_REGION, temp_center$REGION_TYPE, "Telomeric"))))
  } else {
    temp_regions <- outlier_by_chrom %>%
      dplyr::rowwise() %>%
      dplyr::mutate(GENOMIC_REGION = ifelse(CHROM == chrom && MID_BIN > temp_Larm$START_REGION && MID_BIN < temp_Larm$END_REGION, "LEFT_ARM",
                                            ifelse(CHROM == chrom && MID_BIN > temp_Rarm$START_REGION && MID_BIN < temp_Rarm$END_REGION, "RIGHT_ARM",
                                                   ifelse(CHROM == chrom && MID_BIN > temp_center$START_REGION && MID_BIN < temp_center$END_REGION, temp_center$REGION_TYPE, "Telomeric"))))
    count_regions <- dplyr::bind_rows(count_regions, temp_regions)
  }
}
# make non-outlier bins rank to 1+the max outlier rank (for plotting purposes)
count_regions$rank[is.na(count_regions$rank)] <- max(count_regions$rank,na.rm = T)+1

count_regions$STRAIN <- args[1]

# plot
outlier_plot <- ggplot(count_regions)+
  aes(x = MID_BIN/1e6)+
  geom_point( aes(color = rank, y = COUNT)) +
  facet_grid(STRAIN~CHROM, space = "free", scales = "free")+
  scale_color_viridis_c(direction = -1, option = "B")+
  theme_bw(18) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Variant Count")+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray70")) +
  labs(color = "Outlier\nRank")

# save outlier plot
ggsave(outlier_plot, 
       filename = glue::glue("{args[1]}_Outlier_Counts.png"),
       dpi = 300,
       height = 4,
       width = 12)

# save outlier df
write.table(count_regions, file = glue::glue("{args[1]}_Outlier_Counts.tsv"), quote = F, col.names = T, row.names = F, sep = "\t")

