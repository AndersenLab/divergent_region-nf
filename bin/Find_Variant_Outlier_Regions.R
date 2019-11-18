
setwd("~/Dropbox/AndersenLab/LabFolders/Stefan/Manuscripts/CePopGen_Paper/")
# outlier detection in windows
HampelFilter <- function (x, k,t0=3){
  n <- length(x)
  y <- x
  ind <- c()
  L <- 1.4826
  for (i in 1:n) {
    if(i < k+1){
      x0 <- median(c(x[1:(k+i)]))
      S0 <- L * x0
    } else if(i > n-(k+1)){
      x0 <- median(c(x[n:(n-k)]))
      S0 <- L * x0
    } else {
      x0 <- median(x[(i - k):(i + k)])
      S0 <- L * x0
    }
    
    if (x[i] > t0 * S0) {
      y[i] <- x0
      ind <- c(ind, i)
    }
  }
  list(y = y, ind = ind)
}

for(bin_size in c("50k", "100k")){
  # load data
  binned_v <- data.frame()
  for(chrom in c("I","II","III","IV","V","X")){
    temp_v <- data.table::fread(glue::glue("Binned_Varaints/{chrom}.{bin_size}.txt"),col.names = c("CHROM","BIN_START", "VARIANT_COUNT")) %>%
      dplyr::mutate(AVG_BIN_POS = (BIN_START + lead(BIN_START))/2) %>%
      dplyr::mutate(AVG_BIN_POS = ifelse(is.na(AVG_BIN_POS), BIN_START, AVG_BIN_POS)) %>%
      dplyr::mutate(BIN_SIZE = mean(c(lead(BIN_START) - BIN_START),na.rm = T ))
    if(nrow(binned_v) == 0){
      binned_v <- temp_v 
    } else {
      binned_v <- dplyr::bind_rows(binned_v, temp_v)
    }
  }
  
  # apply outlier detection function
  mark_outliers <- data.frame()
  for(chrom in c("I","II","III","IV","V","X")){
    
    temp_outlier <- dplyr::filter(binned_v, CHROM == chrom) %>%
      dplyr::pull(VARIANT_COUNT) %>%
      HampelFilter(., k=20, t0=1)
    
    outlier_df <- dplyr::filter(binned_v, CHROM == chrom) %>%
      dplyr::mutate(adjusted = temp_outlier$y) %>%
      dplyr::mutate(outlier = ifelse(VARIANT_COUNT != adjusted | 
                                       (VARIANT_COUNT < 500 & BIN_SIZE > 90000) |
                                       (VARIANT_COUNT < 200 & BIN_SIZE < 90000)  ,"OUTLIER", "NOT"))
    
    if(nrow(mark_outliers) == 0){
      mark_outliers <- outlier_df 
    } else {
      mark_outliers <- dplyr::bind_rows(mark_outliers, outlier_df)
    }
  }
  
  vc_plot <- ggplot(mark_outliers)+
    aes(x = AVG_BIN_POS/1e6)+
    geom_line(aes(y = VARIANT_COUNT))+
    geom_point( aes(fill = outlier, y = VARIANT_COUNT), shape = 21, size = 3) +
    facet_grid(.~CHROM, space = "free", scales = "free")+
    scale_fill_manual(values = c("cadetblue3","hotpink3"))+
    theme_bw(18) +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "Genomic Position (Mb)", y = "Variant Count")+
    theme(strip.background = element_blank())
  
  ggsave(vc_plot, 
         filename = glue::glue("Binned_Varaints/Outlier_regions_k20_1x_{bin_size}bin.pdf"), 
         height = 6, width = 18)
  
  outlier_bins <- mark_outliers %>%
    dplyr::mutate(outlier_window = ifelse(outlier == "OUTLIER", glue::glue("{BIN_START}-{lead(BIN_START)}"), NA))
  
  write.table(outlier_bins, 
              file = glue::glue("Binned_Varaints/Outlier_bins_k20_1x_{bin_size}.tsv"),
              quote = F, col.names = T, row.names = F, sep = "\t")
}




