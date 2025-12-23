args <- commandArgs(TRUE)
phecode <- as.character(args[1])

##load the package and the function
library(data.table)
source("/net/fantasia/home/borang/MultiAncestry_Coloc/MECOLOC/DataQC.R")
download_UKB_summary <- function(code) {
list=fread("/net/mulan/disk2/luliuu/project4/realdata/Gtex/analysis/Neale_lab_list.txt",header=T)
tmp=list[list$"Phenotype Code"==code,]
cmd=tmp[tmp$Sex=="female",]$"wget command"[1]
system(cmd)
cmd=tmp[tmp$Sex=="male",]$"wget command"[1]
system(cmd)
}

process_UKB_summary <- function(phecode) {
f=list.files(getwd())
GWAS_Male <- fread(f[grep("v3.male.",f)])
GWAS_Female <- fread(f[grep("v3.female.",f)])

GWAS_Male=GWAS_Male[GWAS_Male$low_confidence_variant==FALSE,]
GWAS_Male=na.omit(GWAS_Male)
GWAS_Female=GWAS_Female[GWAS_Female$low_confidence_variant==FALSE,]
GWAS_Female=na.omit(GWAS_Female)

chr_pos_allele<-matrix(unlist(strsplit(GWAS_Male$variant,":")),byrow = T,ncol=4)
GWAS_Male$CHR=chr_pos_allele[,1]
GWAS_Male$POS=as.numeric(chr_pos_allele[,2])
GWAS_Male$REF=chr_pos_allele[,3]
GWAS_Male$ALT=chr_pos_allele[,4]
chr_pos_allele<-matrix(unlist(strsplit(GWAS_Female$variant,":")),byrow = T,ncol=4)
GWAS_Female$CHR=chr_pos_allele[,1]
GWAS_Female$POS=as.numeric(chr_pos_allele[,2])
GWAS_Female$REF=chr_pos_allele[,3]
GWAS_Female$ALT=chr_pos_allele[,4]

GWAS_Male<-GWAS_QC(GWAS_Male,0.001)
GWAS_Female<-GWAS_QC(GWAS_Female,0.001)
GWAS_common<-GWAS_Match(GWAS_Male,GWAS_Female)
chr_pos=paste0(GWAS_common$CHR,"_",GWAS_common$POS)

GWAS_Male<-GWAS_Male[GWAS_Male$CHR_POS %in% chr_pos,]
GWAS_Female<-GWAS_Female[GWAS_Female$CHR_POS %in%chr_pos,]

#GWAS_Male$beta[GWAS_Male$REF!=GWAS_common$REF]=(-1)*GWAS_Male$beta[GWAS_Male$REF!=GWAS_common$REF]
#GWAS_Male$REF=GWAS_common$REF
#GWAS_Male$ALT=GWAS_common$ALT
  
GWAS_Female$beta[GWAS_Female$REF!=GWAS_common$REF]=(-1)*GWAS_Female$beta[GWAS_Female$REF!=GWAS_common$REF]
GWAS_Female$REF=GWAS_common$REF
GWAS_Female$ALT=GWAS_common$ALT

z_score=(GWAS_Female$beta-GWAS_Male$beta)/sqrt((GWAS_Female$se)^2+(GWAS_Male$se)^2)
p=2 * (1 - pnorm(abs(z_score)))
GWAS_common$beta_F=GWAS_Female$beta
GWAS_common$se_F=GWAS_Female$se
GWAS_common$beta_M=GWAS_Male$beta
GWAS_common$se_M=GWAS_Male$se
GWAS_common$sdz=z_score
GWAS_common$sdp=p
GWAS_sd=GWAS_common
save(GWAS_sd,file="GWAS_sd.RDat")
save(GWAS_Male,GWAS_Female,file="GWAS.RDat")

if(sum(GWAS_sd$sdp<5e-8)>0){
sig=GWAS_sd[GWAS_sd$sdp<5e-8,]
write.table(sig,"GWAS_sdsig.txt",quote=FALSE,row.names=FALSE)
}

f=list.files(getwd())
cmd=paste0("rm ",f[grep(".tsv.bgz",f)])
system(cmd[1])
system(cmd[2])
}


setwd("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8")
dir.create(as.character(phecode))
setwd(as.character(phecode))

download_UKB_summary(phecode)

process_UKB_summary(phecode)
