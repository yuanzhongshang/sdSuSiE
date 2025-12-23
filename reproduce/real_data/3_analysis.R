args <- commandArgs(TRUE)
i <- as.integer(args[1])
phecode <- as.character(args[2])
#put this file into the 4079_irnt folder
region=read.table(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/",phecode,"/region.txt"),header=T)

#construct the inputs
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/",phecode))
dir.create("data")
library(data.table)
library("BEDMatrix")
process_geno <- function(geno_data) {
  apply(geno_data, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x - mean(x))
    return(x)
  })
}

chr=region[i,]$chrom
GWAS_Female=fread(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/",phecode,"/GWAS_Female.txt"),header=T)
GWAS_Male=fread(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/",phecode,"/GWAS_Male.txt"),header=T)

ref_M=BEDMatrix(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/MALE/chr_",chr))
ref_F=BEDMatrix(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/FEMALE/chr_",chr))
ref_M_bim<-fread(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/MALE/chr_",chr,".bim"))
ref_F_bim<-fread(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/FEMALE/chr_",chr,".bim"))
ref_M_bim$CHR_POS=paste0(ref_M_bim$V1,"_",ref_M_bim$V4)
ref_F_bim$CHR_POS=paste0(ref_F_bim$V1,"_",ref_F_bim$V4)

GWAS_M=GWAS_Male[GWAS_Male$Region==i,]
GWAS_F=GWAS_Female[GWAS_Female$Region==i,]

GWAS_M=as.data.frame(GWAS_M[,c(3,9,10,11,8)])
colnames(GWAS_M)=c("SNP","Beta","Se","Z","N")

snp=GWAS_M$SNP
idx=which(ref_M_bim$CHR_POS %in% snp)
if(length(idx)>length(GWAS_M$SNP)){
t1=paste0(GWAS_F$CHR_POS,"_",GWAS_F$REF,"_",GWAS_F$ALT)
t2=paste0(GWAS_F$CHR_POS,"_",GWAS_F$ALT,"_",GWAS_F$REF)
t3=paste0(ref_M_bim$CHR_POS,"_",ref_M_bim$V6,"_",ref_M_bim$V5)

idx1=which(t3 %in% t1)
idx2=which(t3 %in% t2)
idx=union(idx1,idx2)
idx=sort(idx)
}

GWAS_F=as.data.frame(GWAS_F[,c(3,9,10,11,8)])
colnames(GWAS_F)=c("SNP","Beta","Se","Z","N")

M=process_geno(ref_M[,idx])
LD_M=Rfast::cora(M)
colnames(LD_M)=GWAS_M$SNP
rownames(LD_M)=GWAS_M$SNP
F=process_geno(ref_F[,idx])
LD_F=Rfast::cora(F)
colnames(LD_F)=GWAS_M$SNP
rownames(LD_F)=GWAS_M$SNP
na_columns1 <- which(colSums(is.na(LD_M)) == (nrow(LD_M)-1))
na_columns2 <- which(colSums(is.na(LD_F)) == (nrow(LD_F)-1))
del=union(na_columns1,na_columns2)
if(length(del)>0){
GWAS_M=GWAS_M[-del,]
GWAS_F=GWAS_F[-del,]
LD_M=LD_M[-del,-del]
LD_F=LD_F[-del,-del]
}
save(GWAS_M,GWAS_F,LD_M,LD_F, file = paste0("data/data_",i,".RData"))

#sdSuSiE
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/",phecode))
dir.create("susieresult")
dir.create("susieresultL_10")
dir.create("susiePIP")
library(sdSuSiE)

summary_list=list()
summary_list[[1]]=GWAS_M
summary_list[[2]]=GWAS_F
names(summary_list)=c("Male","Female")

mat_list=list()
mat_list[[1]]=LD_M
mat_list[[2]]=LD_F
names(mat_list)=c("Male","Female")

our<-sdSuSiE(mat_list,summary_list,L=10)
save(our,file = paste0("susieresultL_10/region_",i,".RData"))

L=10
if (!(length(na.omit(our$ELBO)) != 101) ) {
  L_values <- 9:1
  for (L in L_values) {
    our <- sdSuSiE(mat_list,summary_list,L=L)
    if (length(na.omit(our$ELBO)) != 101) break
  }
}
save(our,file = paste0("susieresult/region_",rep,"_",L,".RData"))
data <- data.frame(
    SNP       = summary_list[[1]]$SNP,
    PIPEither = our$PIP,
    PIPc      = our$PIP_config[,1],
    PIPd      = our$PIP_config[,2],
    PIPb  = our$PIP_config[,3]
  )
write.table(data,paste0("susiePIP/region_",rep,"_",L,".txt"),quote=FALSE)

#stepwise regression
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/",phecode))
dir.create("stepwise")

SNP = GWAS_M$SNP
Z = (sqrt(GWAS_M$N - 1) * GWAS_M$Z + sqrt(GWAS_F$N - 1) * GWAS_F$Z) / sqrt(GWAS_M$N + GWAS_F$N - 1)
N = round(GWAS_M$N + GWAS_F$N)
Beta = Z / sqrt(N - 1)
Se = 1 / sqrt(N - 1)
data1 = data.frame(SNP, Beta, Se, Z, N)

SNP = paste0("S*", SNP)
Z = (-sqrt(GWAS_M$N - 1) * GWAS_M$Z + sqrt(GWAS_F$N - 1) * GWAS_F$Z) / sqrt(GWAS_M$N + GWAS_F$N - 1)
N = round(GWAS_M$N + GWAS_F$N)
Beta = Z / sqrt(N - 1)
Se = 1 / sqrt(N - 1)
data2 = data.frame(SNP, Beta, Se, Z, N)

data = rbind(data1, data2)
data$p = pchisq((data$Beta/data$Se)^2, df = 1, lower.tail = FALSE)

N_M = round(median(GWAS_M$N))
N_F = round(median(GWAS_F$N))
N = N_M + N_F
V1 = (N_M-1)*LD_M + (N_F-1)*LD_F
V2 = (-1)*(N_M-1)*LD_M + (N_F-1)*LD_F
LD = rbind(cbind(V1,V2),cbind(V2,V1)) / (N-1)

res = sdstepwise(data, LD)
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/",phecode,"/stepwise"))
write.table(res,paste0("res_",i,".txt"), quote=FALSE, row.names=FALSE)

#SuSiE-modify
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/",phecode))
dir.create("susie_modify")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/",phecode,"/susie_modify"))
dir.create("susieresult")
dir.create("susiePIP")
library(susieR)

sd=(GWAS_M$Z/sqrt(GWAS_M$N-1)-GWAS_F$Z/sqrt(GWAS_F$N-1))/sqrt(1/(GWAS_M$N-1)+1/(GWAS_F$N-1))
snp=GWAS_F$SNP
N=GWAS_M$N[1]+GWAS_F$N[1]

R=((GWAS_M$N[1]-1)*LD_M+(GWAS_F$N[1]-1)*LD_F)/(N-1)
fit = susie_rss(sd,R,N)
save(fit,file = paste0("susieresult/region_",i,".RData"))
PIP=fit$pip
data=data.frame(snp,PIP)
write.table(data,paste0("susiePIP/region_",i,".txt"),quote=FALSE)

##MESuSiE-modify
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/",phecode))
dir.create("mesusie_modify")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/",phecode,"/mesusie_modify"))
library(MESuSiE)

our<-meSuSie_core(mat_list,summary_list,L=10)
sumsEB1 <- Reduce("+", our$EB1)
beta_EA=sumsEB1[,1]
beta_AA=sumsEB1[,2]

var_EA=0
var_AA=0
for(j in 1:length(our$alpha)){
var_EA=var_EA+our$alpha[[j]][,1]*(our$mu2[[j]][[1]]-our$mu1[[j]][[1]]^2)+our$alpha[[j]][,3]*(our$mu2[[j]][[3]][,1]-our$mu1[[j]][[3]][,1]^2)
var_AA=var_AA+our$alpha[[j]][,2]*(our$mu2[[j]][[2]]-our$mu1[[j]][[2]]^2)+our$alpha[[j]][,3]*(our$mu2[[j]][[3]][,2]-our$mu1[[j]][[3]][,2]^2)
}
z=(beta_EA-beta_AA)/sqrt(var_EA+var_AA)
log_p_value <- pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
sdp_L10 <- 2 * exp(log_p_value)
sdr=data.frame(snp,sdp_L10)
L=10
if (!((length(na.omit(our$ELBO)) != 101) & (our$ELBO[length(na.omit(our$ELBO))]-our$ELBO[length(na.omit(our$ELBO))-1]>0))) {
  L_values <- 9:1
  for (L in L_values) {
    our <- meSuSie_core(mat_list,summary_list,L=L)
    if ((length(na.omit(our$ELBO)) != 101) & (our$ELBO[length(na.omit(our$ELBO))]-our$ELBO[length(na.omit(our$ELBO))-1]>0)) break
  }
}
sumsEB1 <- Reduce("+", our$EB1)
beta_EA=sumsEB1[,1]
beta_AA=sumsEB1[,2]

var_EA=0
var_AA=0
for(j in 1:length(our$alpha)){
var_EA=var_EA+our$alpha[[j]][,1]*(our$mu2[[j]][[1]]-our$mu1[[j]][[1]]^2)+our$alpha[[j]][,3]*(our$mu2[[j]][[3]][,1]-our$mu1[[j]][[3]][,1]^2)
var_AA=var_AA+our$alpha[[j]][,2]*(our$mu2[[j]][[2]]-our$mu1[[j]][[2]]^2)+our$alpha[[j]][,3]*(our$mu2[[j]][[3]][,2]-our$mu1[[j]][[3]][,2]^2)
}
z=(beta_EA-beta_AA)/sqrt(var_EA+var_AA)
log_p_value <- pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
sdp_L <- 2 * exp(log_p_value)
sdr$sdp_L=sdp_L
sdr$L=L
sdr$beta_EA=beta_EA
sdr$beta_AA=beta_AA
sdr$var_EA=var_EA
sdr$var_AA=var_AA

save(sdr, file = paste0("data_",i,".RData"))
