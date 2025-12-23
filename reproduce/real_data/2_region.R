args <- commandArgs(TRUE)
phecode <- as.character(args[1])

setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/",phecode))

##check if exists the sdSNP with p-value<5e-8
if(file.exists("GWAS_sdsig.txt")){
library(data.table)
library(dplyr)
GWAS_sd=fread("GWAS_sdsig.txt",header=T)
GWAS_sd$CHR_POS=paste0(GWAS_sd$CHR,"_",GWAS_sd$POS)
GWAS_sd$CHR[GWAS_sd$CHR=="X"]="23"
signal_SNP=GWAS_sd%>%filter(as.numeric(sdp)<5e-8)%>%select(CHR,POS)

if(dim(signal_SNP)[1]>0){
##load the package and the function
source("/net/fantasia/home/borang/MultiAncestry_Coloc/MECOLOC/DataQC.R")
library(bedr)
library(GenomicRanges)
library(patchwork)

#######Define functions to be used
find_mismatch<-function(dat1,dat2,REF_1,ALT_1,REF_2,ALT_2){
  pos_idx<-which(!((dat1%>%dplyr::select(REF_1)==dat2%>%dplyr::select(REF_2) & dat1%>%dplyr::select(ALT_1)==dat2%>%dplyr::select(ALT_2))|(dat1%>%dplyr::select(REF_1)==dat2%>%dplyr::select(ALT_2)&dat1%>%dplyr::select(ALT_1)==dat2%>%dplyr::select(REF_2))))
  return(pos_idx)
}

sumstat_ref<-function(sumstat,ref){
	common_SNP<-intersect(sumstat$CHR_POS,ref$CHR_POS)
	sumstat_sub<-sumstat%>%filter(CHR_POS%in%common_SNP)%>%arrange(match(CHR_POS, common_SNP))
	ref_sub<-ref%>%filter(CHR_POS%in%common_SNP)%>%arrange(match(CHR_POS, common_SNP))
	pos_idx<-find_mismatch(sumstat_sub,ref_sub,"REF","ALT","V6","V5")
	if(length(pos_idx)!=0)common_SNP<-common_SNP[-pos_idx]
	return(common_SNP)
}
range_to_region<-function(dat,common_SNP){
  Reduce(rbind,lapply(1:nrow(dat),function(region){
    chr<-dat$chrom[region]
    start<-dat$start[region]
    stop<-dat$stop[region]
    region_dat<-common_SNP%>%filter(CHR==chr,POS>start,POS<stop)
    region_dat$region = region
    cat(region)
    return(region_dat)
  }))
}

bim_ukb=fread("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/bim_ukb.txt")
bim_ukb=as_tibble(bim_ukb)

load("GWAS.RDat")
GWAS_Male$CHR[GWAS_Male$CHR=="X"]="23"
GWAS_Male$CHR_POS[GWAS_Male$CHR=="23"]=gsub("X","23",GWAS_Male$CHR_POS[GWAS_Male$CHR=="23"])
GWAS_Female$CHR[GWAS_Female$CHR=="X"]="23"
GWAS_Female$CHR_POS[GWAS_Female$CHR=="23"]=gsub("X","23",GWAS_Female$CHR_POS[GWAS_Female$CHR=="23"])
bim_ukb=bim_ukb[!duplicated(bim_ukb$CHR_POS),]
common_SNP_UKB<-sumstat_ref(GWAS_Male,bim_ukb)
GWAS_Male<-GWAS_Male[match(common_SNP_UKB,GWAS_Male$CHR_POS),]
GWAS_Female<-GWAS_Female[match(common_SNP_UKB,GWAS_Female$CHR_POS),]
bim_ukb<-bim_ukb[match(common_SNP_UKB,bim_ukb$CHR_POS),]

signal_SNP=signal_SNP[paste0(signal_SNP$CHR,"_",signal_SNP$POS) %in% common_SNP_UKB,]

  signal_SNP_update<-c()
    for(chr in 1:22){
      signal_SNP_subset<-signal_SNP[signal_SNP$CHR==chr,]
	  bim_ukb_subset<-bim_ukb[bim_ukb$V1==chr,]
      signal_SNP_subset$start<-signal_SNP_subset$POS-5e5
      signal_SNP_subset$start<-ifelse(signal_SNP_subset$start<min(bim_ukb_subset$V4),min(bim_ukb_subset$V4),signal_SNP_subset$start)      
      signal_SNP_subset$end<-signal_SNP_subset$POS+5e5
      signal_SNP_subset$end<-ifelse(signal_SNP_subset$end>max(bim_ukb_subset$V4),max(bim_ukb_subset$V4),signal_SNP_subset$end)
      signal_SNP_update<-rbind(signal_SNP_update,signal_SNP_subset)
    }
	chr=23
	  signal_SNP_subset<-signal_SNP[signal_SNP$CHR==chr,]
	  bim_ukb_subset<-bim_ukb[bim_ukb$V1==chr,]
      signal_SNP_subset$start<-signal_SNP_subset$POS-5e5
      signal_SNP_subset$start<-ifelse(signal_SNP_subset$start<min(bim_ukb_subset$V4),min(bim_ukb_subset$V4),signal_SNP_subset$start)      
      signal_SNP_subset$end<-signal_SNP_subset$POS+5e5
      signal_SNP_subset$end<-ifelse(signal_SNP_subset$end>max(bim_ukb_subset$V4),max(bim_ukb_subset$V4),signal_SNP_subset$end)
      signal_SNP_update<-rbind(signal_SNP_update,signal_SNP_subset)
    signal_SNP<-signal_SNP_update
    chr_region<-paste0("chr",signal_SNP$CHR,":",signal_SNP$start,"-",signal_SNP$end)
    is.sorted <- is.sorted.region(chr_region)
    chr_region.sort.natural <- bedr.sort.region(chr_region, method = "natural")
    is.merged <- is.merged.region(chr_region.sort.natural)
    chr_region.sort.natural.merge <- bedr.merge.region(chr_region.sort.natural)
    candidate_region<-Reduce(rbind,lapply(1:length(chr_region.sort.natural.merge),function(region){
      chr<-as.numeric(gsub("chr","",gsub(":.*","",chr_region.sort.natural.merge[region])))
      start<-as.numeric(gsub("-.*","",gsub(".*:","",chr_region.sort.natural.merge[region])))
      end<-as.numeric(gsub(".*-","",gsub(".*:","",chr_region.sort.natural.merge[region])))
      return(c(chr,start,end))
    }))
	if(length(candidate_region)==3){
	candidate_region=t(matrix(candidate_region))
}
    candidate_region<-data.frame(chrom = candidate_region[,1],start = candidate_region[,2],stop = candidate_region[,3])
	candidate_region=candidate_region[order(candidate_region$chrom),]
write.table(candidate_region,"region.txt",quote=FALSE, row.names=FALSE)

#Select all SNPs that are in the candidate region
candidate_region_SNP<-range_to_region(candidate_region,GWAS_Male)
GWAS_Male_pro<-GWAS_Male[match(candidate_region_SNP$CHR_POS,GWAS_Male$CHR_POS),]
GWAS_Male_pro<-GWAS_Male_pro%>%mutate(Region = candidate_region_SNP$region)
	
GWAS_Female_pro<-GWAS_Female[match(candidate_region_SNP$CHR_POS,GWAS_Female$CHR_POS),]
GWAS_Female_pro<-GWAS_Female_pro%>%mutate(Region = candidate_region_SNP$region)

#Write out the processed GWAS sumstat in cnadidate region
GWAS_Female_pro<-dplyr::select(GWAS_Female_pro,c(CHR,POS,CHR_POS,REF,ALT,minor_allele,minor_AF,n_complete_samples,beta,se,tstat,pval,Region))
GWAS_Male_pro<-dplyr::select(GWAS_Male_pro,c(CHR,POS,CHR_POS,REF,ALT,minor_allele,minor_AF,n_complete_samples,beta,se,tstat,pval,Region))
bim_ref=bim_ukb[bim_ukb$CHR_POS %in% GWAS_Female_pro$CHR_POS,]

#covert the mismatched allele
idx=which(GWAS_Female_pro$REF!=bim_ref$V6)
GWAS_Male_pro$beta[idx]<-(-1)*GWAS_Male_pro$beta[idx]
GWAS_Male_pro$tstat[idx]<-(-1)*GWAS_Male_pro$tstat[idx]
GWAS_Male_pro$REF[idx]<-bim_ref$V6[idx]
GWAS_Male_pro$ALT[idx]<-bim_ref$V5[idx]
GWAS_Female_pro$beta[idx]<-(-1)*GWAS_Female_pro$beta[idx]
GWAS_Female_pro$tstat[idx]<-(-1)*GWAS_Female_pro$tstat[idx]
GWAS_Female_pro$REF[idx]<-bim_ref$V6[idx]
GWAS_Female_pro$ALT[idx]<-bim_ref$V5[idx]

write.table(GWAS_Female_pro,"GWAS_Female.txt",quote=FALSE,row.names=FALSE)
write.table(GWAS_Male_pro,"GWAS_Male.txt",quote=FALSE,row.names=FALSE)
}
}
