args <- commandArgs(TRUE)
rep <- as.integer(args[1])

regionanaly=read.table("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/simulation_region.txt")

chr=regionanaly[rep,1]
start=regionanaly[rep,2]
stop=regionanaly[rep,3]

set.seed(1217-rep)

library(data.table)
library(sdSuSiE)
library(mvtnorm)
library("BEDMatrix")

#load the pre-calculated extracted genotype and LD matrix across sexes
load(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/geno/region_",rep,".RData"))
load(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LD/region_",rep,".RData"))

process_geno <- function(geno_data) {
  apply(geno_data, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x - mean(x))
    return(x)
  })
}
ref_M_bim<-fread(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/MALE/chr_",chr,".bim"))
ref_M_bim$CHR_POS=paste0(ref_M_bim$V1,"_",ref_M_bim$V4)
snp <- sub("_.*", "", colnames(ref_M_spe))
target=ref_M_bim[match(snp,ref_M_bim$V2),]
target$CHR_POS_A1_A2=paste0(target$CHR_POS,"_",target$V5,"_",target$V6)
target$CHR_POS_A2_A1=paste0(target$CHR_POS,"_",target$V6,"_",target$V5)
if(chr<23){
d1000G=BEDMatrix(paste0("/net/mulan/boran/geno/1000G/plink/EUR_filter/chr",chr))
d1000G_bim<-fread(paste0("/net/mulan/boran/geno/1000G/plink/EUR_filter/chr",chr,".bim"))
d1000G_fam<-fread(paste0("/net/mulan/boran/geno/1000G/plink/EUR_filter/chr",chr,".fam"))
}else{
d1000G=BEDMatrix("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/1000Gdata/chrX/chrX")
d1000G_bim<-fread("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/1000Gdata/chrX/chrX.bim")
d1000G_fam<-fread("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/1000Gdata/chrX/chrX.fam")
d1000G_bim$V1=23
}
d1000G_bim$CHR_POS_A1_A2=paste0(d1000G_bim$V1,"_",d1000G_bim$V4,"_",d1000G_bim$V5,"_",d1000G_bim$V6)
d1000G_bim$CHR_POS_A2_A1=paste0(d1000G_bim$V1,"_",d1000G_bim$V4,"_",d1000G_bim$V6,"_",d1000G_bim$V5)
sex=fread("/net/mulan/boran/geno/1000G/vcf/integrated_call_samples_v3.20130502.ALL.panel")
sex=sex[match(d1000G_fam$V2,sex$sample),]
m_idx=which(sex$gender=="male")
f_idx=which(sex$gender=="female")

sel=which(target$CHR_POS_A1_A2 %in% d1000G_bim$CHR_POS_A1_A2 | target$CHR_POS_A2_A1 %in% d1000G_bim$CHR_POS_A2_A1)
LD_M=LD_M[sel,sel]
LD_F=LD_F[sel,sel]
ref_M_spe=ref_M_spe[,sel]
ref_F_spe=ref_F_spe[,sel]
target=target[sel,]
sel=which(d1000G_bim$CHR_POS_A1_A2 %in% target$CHR_POS_A1_A2 | d1000G_bim$CHR_POS_A1_A2 %in% target$CHR_POS_A2_A1)
d1000G_bim=d1000G_bim[sel,]
d1000G=d1000G[,sel]
if(sum(d1000G_bim$V5!=target$V5)){
idx=which(d1000G_bim$V5!=target$V5)
d1000G <- d1000G %>%
  mutate(across(all_of(names(d1000G)[idx]), ~ ifelse(. == 0, 2, ifelse(. == 2, 0, .))))
}

invariant_cols_m <- apply(d1000G[m_idx,], 2, function(col) length(unique(col)) == 1)
invariant_cols_f <- apply(d1000G[f_idx,], 2, function(col) length(unique(col)) == 1)
invariant_cols=sort(union(which(invariant_cols_m==1),which(invariant_cols_f==1)))
LD_M=LD_M[-invariant_cols,-invariant_cols]
LD_F=LD_F[-invariant_cols,-invariant_cols]
ref_M_spe=ref_M_spe[,-invariant_cols]
ref_F_spe=ref_F_spe[,-invariant_cols]
target=target[-invariant_cols,]
d1000G=d1000G[,-invariant_cols]

LD_M_1000G=Rfast::cora(process_geno(d1000G[m_idx,]))
LD_F_1000G=Rfast::cora(process_geno(d1000G[f_idx,]))

#set the parameter
h_c=5e-4
h_d=5e-4
n_causal=1
ref_N_M=50000
ref_N_F=50000
N_M=150000
N_F=150000
rho=0.5

select_cd=sort(sample(c(1:dim(LD_M)[1]),n_causal))
beta_c=rep(0,dim(LD_M)[1])
beta_d=rep(0,dim(LD_M)[1])

beta_shared<-rmvnorm(length(select_cd),mean=rep(0,2),sigma=matrix(c(h_c,rho*sqrt(h_d*h_c),rho*sqrt(h_d*h_c),h_d)/n_causal,ncol=2,nrow=2))
beta_c[select_cd]=beta_shared[,1]
beta_d[select_cd]<-beta_shared[,2]

M_marginal<-LD_M%*%(beta_c-beta_d)
F_marginal<-LD_F%*%(beta_c+beta_d)

y_null_M<-rnorm(ref_N_M,0,sqrt(1-h_c-h_d))
err_beta_M<-t(ref_M_spe)%*%y_null_M/ref_N_M
y_null_F<-rnorm(ref_N_F,0,sqrt(1-h_c-h_d))
err_beta_F<-t(ref_F_spe)%*%y_null_F/ref_N_F

err_beta_M_scale<-sqrt(ref_N_M/N_M)*err_beta_M
err_beta_F_scale<-sqrt(ref_N_F/N_F)*err_beta_F
      
z_M<-(M_marginal+err_beta_M_scale)*sqrt(N_M)
z_F<-(F_marginal+err_beta_F_scale)*sqrt(N_F)  

setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LD1000G/"))
dir.create("causal")
save(select_cd,file = paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LD1000G/causal/data_",rep,".RData"))

snp=target$V2

##sdsusie
Z=z_M
Beta=Z/sqrt(N_M-1)
Se=1/sqrt(N_M-1)
N=N_M
data1=data.frame(snp,Beta,Se,Z,N)
colnames(data1)[1]="SNP"
Z=z_F
Beta=Z/sqrt(N_F-1)
Se=1/sqrt(N_F-1)
N=N_F
data2=data.frame(snp,Beta,Se,Z,N)
colnames(data2)[1]="SNP"

summary_list=list()
summary_list[[1]]=data1
summary_list[[2]]=data2
names(summary_list)=c("Male","Female")

mat_list=list()
mat_list[[1]]=LD_M
mat_list[[2]]=LD_F
names(mat_list)=c("Male","Female")

dir.create("gold")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LD1000G/gold"))
dir.create("susieresult")
dir.create("susieresultL_10")
dir.create("susiePIP")

our<-sdSuSiE(mat_list,summary_list,L=10)
save(our,file = paste0("susieresultL_10/region_",rep,".RData"))

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
    PIPShare  = our$PIP_config[,3]
  )
write.table(data,paste0("susiePIP/region_",rep,"_",L,".txt"),quote=FALSE)

#1000G
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LD1000G/"))
dir.create("1000G")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LD1000G/1000G"))

mat_list=list()
mat_list[[1]]=LD_M_1000G
mat_list[[2]]=LD_F_1000G
names(mat_list)=c("Male","Female")

dir.create("susie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LD1000G/1000G/susie"))
dir.create("susieresult")
dir.create("susieresultL_10")
dir.create("susiePIP")

our<-sdSuSiE(mat_list,summary_list,L=10)
save(our,file = paste0("susieresultL_10/region_",rep,".RData"))

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
    PIPShare  = our$PIP_config[,3]
  )
write.table(data,paste0("susiePIP/region_",rep,"_",L,".txt"),quote=FALSE)

