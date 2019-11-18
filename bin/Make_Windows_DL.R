library(tidyverse)
library(zoo)
library(rlang)

options(scipen=999)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

# args
# 1 - processed outlier file
# 2 - fraction of window to be overlapped by SV (--sv_fraction)
# 3 - minimum variant count to be considered diverged (determined by --percent_outlier in previous process)
# 4 - number of SVs for SV coverage mask (--nsv_svcov)
# 5 - coverage for SV coverage mask (--cov_svcov)
# example
# args <- c("CB4856_Processed_Outliers.tsv", "0.5", "22", "5", "5")

args <- commandArgs(TRUE)

# leads <- function(var, n=as.numeric(args[2])){
#   var <- enquo(var)
#   
#   indices <- seq_len(n)
#   map( indices, ~quo(dplyr::lead(!!var, !!.x)) ) %>% 
#     set_names(sprintf("lag_%s_%02d", quo_text(var), indices))
#   
# }
# 
# lags <- function(var, n=as.numeric(args[2])){
#   var <- enquo(var)
#   
#   indices <- seq_len(n)
#   map( indices, ~quo(dplyr::lag(!!var, !!.x)) ) %>% 
#     set_names(sprintf("lag_%s_%02d", quo_text(var), indices))
#   
# }
# read data
outliers <- readr::read_tsv(args[1])

# makin windows
m_windows <- outliers %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN, .keep_all = T) %>%
  dplyr::group_by(CHROM) %>%
  dplyr::mutate(window_mask = ifelse(cov_direction %in% c("Down"), "Masked_Low_Coverage", "Pass")) %>% # coverage mask
  dplyr::mutate(window_mask = ifelse(COUNT > as.numeric(args[3]), ifelse(window_mask == "Pass", "Masked_Count", paste(window_mask, "Masked_Count", sep="&")), window_mask)) %>% # variant count mask
  dplyr::mutate(window_mask = ifelse(fraction_SV_bases > as.numeric(args[2]), ifelse(window_mask == "Pass", "Masked_SV", paste(window_mask, "Masked_SV", sep="&")), window_mask)) %>% # High SV fraction
  dplyr::mutate(window_mask = ifelse(fraction_SV_bases > 0 & COUNT > as.numeric(args[3])*(1-fraction_SV_bases), ifelse(window_mask == "Pass", "Masked_SV_count", paste(window_mask, "Masked_SV_count", sep="&")), window_mask)) # SVs and adjusted variant count

# Counts
m_windows_pr <- m_windows %>%
  # If regions flanking are both masked, then mask central region
  dplyr::mutate(window_mask = ifelse((grepl("Masked", dplyr::lag(window_mask)) & grepl("Masked", dplyr::lead(window_mask))) & !grepl("Masked", window_mask), "Masked_Two_Flank", window_mask))

masked_plot <- ggplot(m_windows_pr)+
  facet_grid(CHROM~., space = "free", scales = "free")+
  geom_point(aes(x = MID_BIN/1e6, y = COUNT, color = window_mask)) +
  geom_segment(aes(x = START_BIN/1e6, xend = END_BIN/1e6, y = -5, yend = -5, color = window_mask), size = 3,
               data = m_windows_pr %>% dplyr::filter(grepl("Mask", window_mask)))+
  theme_bw(18) +
  scale_color_viridis_d(direction = -1)+
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Variant Count")+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray80"))

strain <- unique(m_windows$STRAIN)

ggsave(masked_plot, filename = glue::glue("{strain}_masked_plot.pdf"), height = 10, width = 24)

write.table(m_windows_pr, file = glue::glue("{strain}_Mask_DF.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")

loose_mask <- m_windows_pr %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("Mask", window_mask)) %>%
  dplyr::select(CHROM, START_BIN, END_BIN) %>%
  tidyr::unite(temp1, START_BIN,END_BIN, sep = "-") %>%
  tidyr::unite(REGION, CHROM, temp1, sep = ":")

write.table(loose_mask, file = glue::glue("{strain}_Mask_Regions.tsv"), col.names = F, row.names = F, quote = F, sep = "\t")

