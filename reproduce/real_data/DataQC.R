
require(dplyr)
#' Quality Control for GWAS Data
#'
#' Performs quality control (QC) on GWAS summary statistics to remove SNPs based on certain criteria.
#' Note that the QC control steps are based on the summary statistics have already passed genotype call rate, INFO, and HWE
#'
#' @param dat A data frame containing GWAS summary statistics. Must have columns: CHR, POS, REF, ALT, and MAF.
#' @param MAF_threshold A numeric threshold for minor allele frequencies (MAF) below which SNPs will be removed. Default is 0.01.
#' @param MHC_rm A logical flag indicating whether to remove SNPs located within the MHC region on chromosome 6. Default is TRUE.
#' 
#' @return A data frame after applying QC filters.
#' @export
GWAS_QC <- function(dat, MAF_threshold = 0.01, MHC_rm = TRUE) {
  
  dat<-dat%>%mutate(CHR_POS = paste0(CHR,"_",POS))%>%group_by(CHR_POS) %>%filter(n() == 1) %>% ungroup()
  
  
  # Strand ambiguous
  pos_idx <- which((dat$REF == "G" & dat$ALT == "C") |
                     (dat$REF == "C" & dat$ALT == "G") |
                     (dat$REF == "T" & dat$ALT == "A") |
                     (dat$REF == "A" & dat$ALT == "T"))
  if (length(pos_idx) != 0) dat <- dat[-pos_idx,]
  
  # Multi-allelic
  multi_allelic_index <- dat %>% 
    summarise(N_allele = nchar(REF) + nchar(ALT)) %>% 
    summarise(multi_allelic_index = which(N_allele > 2)) %>% 
    pull(multi_allelic_index)
  if (length(multi_allelic_index) != 0) dat <- dat[-multi_allelic_index,]
  
  #dat<-dat%>%mutate(CHR_POS = paste0(CHR,"_",POS))%>%distinct(CHR_POS,.keep_all=TRUE)%>%select(-CHR_POS)
   
  # Optionally, remove SNPs in the MHC region
  if (MHC_rm) {
    dat <- dat %>% 
      filter(!(CHR == 6 & POS > 25e6 & POS < 34e6))
  }
 
  # MAF
  if("MAF"%in%colnames(dat)){
  dat <- dat %>% 
    filter(MAF > MAF_threshold, MAF < 1 - MAF_threshold)
  }
  
  return(dat)
}


###When there are more than two GWAS to be matched up, we can do something like match the first to all the rest and then find the final set of 

#' Match and Filter GWAS Datasets Based on SNP Identifiers and Alleles
#'
#' This function matches two GWAS datasets based on chromosome and position of SNPs. It then filters the datasets 
#' to retain only common SNPs with matching reference and alternative alleles between the datasets.
#'
#' @param dat1 A data frame containing GWAS data with at least the following columns: CHR, POS, REF, and ALT.
#' @param dat2 A data frame containing GWAS data with at least the following columns: CHR, POS, REF, and ALT.
#'
#' @return A data frame containing the matched and filtered rows from `dat1`.
#'
#' @examples
#' \dontrun{
#'   # Assuming df1 and df2 are your GWAS datasets
#'   matched_data <- GWAS_Match(df1, df2)
#' }
#' 
#' @export
GWAS_Match <- function(dat1, dat2) {
  # Create a unique identifier for each SNP based on CHR and POS
  dat1 <- dat1 %>% mutate(CHR_POS = paste0(CHR, "_", POS))
  dat2 <- dat2 %>% mutate(CHR_POS = paste0(CHR, "_", POS))
  
  # Identify common SNPs based on CHR_POS
  common_SNP <- intersect(dat1$CHR_POS, dat2$CHR_POS)
  
  # Subset each dataset to only include common SNPs
  dat1 <- dat1 %>% 
    filter(CHR_POS %in% common_SNP) %>%
    arrange(match(CHR_POS, common_SNP))
  
  dat2 <- dat2 %>% 
    filter(CHR_POS %in% common_SNP) %>%
    arrange(match(CHR_POS, common_SNP))
  
  # Check if the REF and ALT alleles match between the two datasets
  dat1 <- dat1 %>% mutate(REF_ALT = paste0(REF, ALT), ALT_REF = paste0(ALT, REF))
  dat2 <- dat2 %>% mutate(REF_ALT = paste0(REF, ALT), ALT_REF = paste0(ALT, REF))
  
  # Identify mismatches between the two datasets
  mismatch_idx <- which(!(dat1$REF_ALT == dat2$REF_ALT | dat1$REF_ALT == dat2$ALT_REF))
  
  # Remove rows with mismatches
  if (length(mismatch_idx) > 0) {
    dat1 <- dat1[-mismatch_idx, ]
    dat2 <- dat2[-mismatch_idx, ]
  }
  dat1<-dat1%>%select(CHR,POS,REF,ALT)
  # Return the cleaned dat1
  return(dat1)
}


range_to_region<-function(candidate_region,GWAS){
  Reduce(rbind,lapply(1:nrow(candidate_region),function(region){
    chr<-candidate_region$chrom[region]
    start<-candidate_region$start[region]
    stop<-candidate_region$stop[region]
    region_candidate_region<-GWAS%>%filter(CHR==chr,POS>=start,POS<=stop)
    region_candidate_region$region = region
    cat(region)
    return(region_candidate_region)
  }))
}