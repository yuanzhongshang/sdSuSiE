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

#load the pre-calculated extracted genotype across sexes
load(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/geno/region_",rep,".RData"))
load(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LD/region_",rep,".RData"))

process_geno <- function(geno_data) {
  apply(geno_data, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x - mean(x))
    return(x)
  })
}

h_c=5e-4
h_d=5e-4
n_causal=1
ref_N_M=50000
ref_N_F=50000
N_M=150000
N_F=150000
rho=0.5

rm(select_cd)
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

setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/"))
dir.create("causal")
save(select_cd,file = paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/causal/data_",rep,".RData"))

sd=(z_M/sqrt(N_M-1)-z_F/sqrt(N_F-1))/sqrt(1/(N_M-1)+1/(N_F-1))
log_p_value <- pnorm(abs(sd), lower.tail = FALSE, log.p = TRUE)
sdp <- 2 * exp(log_p_value)
snp=sub("_[^_]*$", "", rownames(sdp))
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

#switched
dir.create("switched")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/switched"))
mat_list=list()
mat_list[[1]]=LD_F
mat_list[[2]]=LD_M
names(mat_list)=c("Male","Female")
dir.create("susie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/switched/susie"))
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

#combined
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/"))
dir.create("combined")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/combined"))

ref_M=BEDMatrix(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/MALE/chr_",chr))
ref_F=BEDMatrix(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/FEMALE/chr_",chr))
ref_M_bim<-fread(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/MALE/chr_",chr,".bim"))
ref_F_bim<-fread(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/FEMALE/chr_",chr,".bim"))
ref_M_bim$CHR_POS=paste0(ref_M_bim$V1,"_",ref_M_bim$V4)
ref_F_bim$CHR_POS=paste0(ref_F_bim$V1,"_",ref_F_bim$V4)

idx=which(ref_F_bim$V4>=start & ref_F_bim$V4<stop)
ref_M_bim=ref_M_bim[idx,]
ref_F_bim=ref_F_bim[idx,]
ref=rbind(as.matrix(ref_M[,idx]),as.matrix(ref_F[,idx]))

ref=process_geno(ref)
LD=Rfast::cora(ref)
mat_list=list()
mat_list[[1]]=LD
mat_list[[2]]=LD
names(mat_list)=c("Male","Female")
dir.create("susie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/combined/susie"))
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

#10% sub-sampled
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/"))
dir.create("01sub")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/01sub"))

LD_M_10=Rfast::cora(process_geno(ref_M[1:5000,idx]))
LD_F_10=Rfast::cora(process_geno(ref_F[1:5000,idx]))
mat_list=list()
mat_list[[1]]=LD_M_10
mat_list[[2]]=LD_F_10
names(mat_list)=c("Male","Female")
dir.create("susie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/01sub/susie"))
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

#0.5% sub-sampled
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/"))
dir.create("0005sub")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/0005sub"))

LD_M_05=Rfast::cora(process_geno(ref_M[1:250,idx]))
LD_F_05=Rfast::cora(process_geno(ref_F[1:250,idx]))

mat_list=list()
mat_list[[1]]=LD_M_05
mat_list[[2]]=LD_F_05
names(mat_list)=c("Male","Female")
dir.create("susie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LDc/0005sub/susie"))
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

