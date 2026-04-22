args <- commandArgs(TRUE)
rep <- as.integer(args[1])

regionanaly=read.table("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/region.txt")

chr=regionanaly[rep,1]
start=regionanaly[rep,2]
stop=regionanaly[rep,3]

library(data.table)
library(sdSuSiE)
library(mvtnorm)

#load the pre-calculated extracted genotype across sexes and LD matrix across sexes
load(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/LD/region_",rep,".RData"))
load(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/geno/region_",rep,".RData"))

#high LD intensity
set.seed(1217-rep)
select_cd=which.max(colSums(LD_M^2))
setwd("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore")
dir.create("high")

#set the parameter
h_c=c(1e-4,5e-4,1e-3)
h_d=5e-4
n_causal=1
ref_N_M=50000
ref_N_F=50000
N_M=150000
N_F=150000
rho=0.5

for(i in 1:3){
setwd("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high")
dir.create(as.character(h_c[i]))
h_d=h_c[i]
beta_c=rep(0,dim(LD_M)[1])
beta_d=rep(0,dim(LD_M)[1])

beta_shared<-rmvnorm(length(select_cd),mean=rep(0,2),sigma=matrix(c(h_c[i],rho*sqrt(h_d*h_c[i]),rho*sqrt(h_d*h_c[i]),h_d)/n_causal,ncol=2,nrow=2))
beta_c[select_cd]=beta_shared[,1]
beta_d[select_cd]<-beta_shared[,2]

M_marginal<-LD_M%*%(beta_c-beta_d)
F_marginal<-LD_F%*%(beta_c+beta_d)

y_null_M<-rnorm(ref_N_M,0,sqrt(1-h_c[i]-h_d))
err_beta_M<-t(ref_M_spe)%*%y_null_M/ref_N_M
y_null_F<-rnorm(ref_N_F,0,sqrt(1-h_c[i]-h_d))
err_beta_F<-t(ref_F_spe)%*%y_null_F/ref_N_F

err_beta_M_scale<-sqrt(ref_N_M/N_M)*err_beta_M
err_beta_F_scale<-sqrt(ref_N_F/N_F)*err_beta_F
      
z_M<-(M_marginal+err_beta_M_scale)*sqrt(N_M)
z_F<-(F_marginal+err_beta_F_scale)*sqrt(N_F)  

setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i]))
dir.create("causal")
save(select_cd,file = paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i],"/causal/data_",rep,".RData"))

#univariate sex-dimorphic analysis
dir.create("sd")
sd=(z_M/sqrt(N_M-1)-z_F/sqrt(N_F-1))/sqrt(1/(N_M-1)+1/(N_F-1))
log_p_value <- pnorm(abs(sd), lower.tail = FALSE, log.p = TRUE)
sdp <- 2 * exp(log_p_value)
snp=sub("_[^_]*$", "", rownames(sdp))
sdr=data.frame(snp,sdp)
save(sdr, file = paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i],"/sd/data_",rep,".RData"))

#stepwise regression
ref_M_bim<-fread(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/MALE/chr_",chr,".bim"))
data=ref_M_bim[match(snp,ref_M_bim$V2),][,c(2,5,6)]
colnames(data)=c("SNP","A1","A2")
N=N_M+N_F
data1=data
data1$b=(sqrt(N_M-1)*z_M+sqrt(N_F-1)*z_F)/(N-1)
data1$se=1/sqrt(N-1)
log_p_value <- pnorm(abs(data1$b/data1$se), lower.tail = FALSE, log.p = TRUE)
data1$p=2 * exp(log_p_value)
data1$N=N
data2=data
data2$b=(-sqrt(N_M-1)*z_M+sqrt(N_F-1)*z_F)/(N-1)
data2$se=1/sqrt(N-1)
log_p_value <- pnorm(abs(data2$b/data2$se), lower.tail = FALSE, log.p = TRUE)
data2$p=2 * exp(log_p_value)
data2$N=N
data2$SNP=paste0("S*",data2$SNP)
data=rbind(data1,data2)
N=N_M+N_F
V1=(N_M-1)*LD_M+(N_F-1)*LD_F
V2=(-1)*(N_M-1)*LD_M+(N_F-1)*LD_F
LD=rbind(cbind(V1,V2),cbind(V2,V1))/(N-1)
res=sdstepwise(data, LD, p_cutoff=5e-8, collinear=0.9,ncore=10)
dir.create("stepwise")
write.table(res,paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i],"/stepwise/res_",rep,".txt"), quote=FALSE, row.names=FALSE)

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

dir.create("susie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i],"/susie"))
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

#SuSiE-modify
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i]))
dir.create("susieadd")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i],"/susieadd"))
dir.create("susieresult")
dir.create("susiePIP")

sd=(z_M/sqrt(N_M-1)-z_F/sqrt(N_F-1))/sqrt(1/(N_M-1)+1/(N_F-1))
snp=sub("_[^_]*$", "", rownames(sd))
N=N_M+N_F
library(susieR)
R=((N_M-1)*LD_M+(N_F-1)*LD_F)/(N-1)
fit = susie_rss(c(sd),R,N)
save(fit,file = paste0("susieresult/region_",rep,".RData"))
PIP=fit$pip
data=data.frame(snp,PIP)
write.table(data,paste0("susiePIP/region_",rep,".txt"),quote=FALSE)

#MESuSiE-modify
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i]))
dir.create("mesusie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i],"/mesusie"))
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
sdp <- 2 * exp(log_p_value)
sdr=data.frame(snp,sdp)
save(sdr, file = paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i],"/mesusie/data_",rep,".RData"))

#SharePro_gxe
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i]))
dir.create("sharepro05")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i],"/sharepro05"))
SNP=sub("_[^_]*$", "", rownames(z_M))
N=N_M
Z=z_M
BETA=Z/sqrt(N_M-1)
SE=1/sqrt(N_M-1)
data1=data.frame(SNP,N,BETA,SE)

Z=z_F
BETA=Z/sqrt(N_F-1)
SE=1/sqrt(N_F-1)
N=N_F
data2=data.frame(SNP,N,BETA,SE)

write.table(data1,paste0("data1_",rep,".txt"),quote=F,row.names=F)
write.table(data2,paste0("data2_",rep,".txt"),quote=F,row.names=F)

write.table((LD_M+LD_F)/2,paste0("LD_",rep,".ld"),quote=F,row.names=F,col.names=F)

system(paste0("python /net/mulan/home/luliuu/SharePro_gxe/src/SharePro/sharepro_gxe.py --z data1_",rep,".txt data2_",rep,".txt --pthres 0.5 --ld LD_",rep,".ld --save res_",rep))
system(paste0("rm data1_",rep,".txt"))
system(paste0("rm data2_",rep,".txt"))
system(paste0("rm LD_",rep,".ld"))
}

#moderate LD intensity
set.seed(1217-rep)
ld_score=colSums(LD_M^2)
select_cd=order(ld_score)[ceiling(length(ld_score)/2)]
setwd("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore")
dir.create("moderate")

#set the parameter
h_c=c(1e-4,5e-4,1e-3)
h_d=5e-4
n_causal=1
ref_N_M=50000
ref_N_F=50000
N_M=150000
N_F=150000
rho=0.5

for(i in 1:3){
setwd("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/moderate")
dir.create(as.character(h_c[i]))
h_d=h_c[i]
beta_c=rep(0,dim(LD_M)[1])
beta_d=rep(0,dim(LD_M)[1])

beta_shared<-rmvnorm(length(select_cd),mean=rep(0,2),sigma=matrix(c(h_c[i],rho*sqrt(h_d*h_c[i]),rho*sqrt(h_d*h_c[i]),h_d)/n_causal,ncol=2,nrow=2))
beta_c[select_cd]=beta_shared[,1]
beta_d[select_cd]<-beta_shared[,2]

M_marginal<-LD_M%*%(beta_c-beta_d)
F_marginal<-LD_F%*%(beta_c+beta_d)

y_null_M<-rnorm(ref_N_M,0,sqrt(1-h_c[i]-h_d))
err_beta_M<-t(ref_M_spe)%*%y_null_M/ref_N_M
y_null_F<-rnorm(ref_N_F,0,sqrt(1-h_c[i]-h_d))
err_beta_F<-t(ref_F_spe)%*%y_null_F/ref_N_F

err_beta_M_scale<-sqrt(ref_N_M/N_M)*err_beta_M
err_beta_F_scale<-sqrt(ref_N_F/N_F)*err_beta_F
      
z_M<-(M_marginal+err_beta_M_scale)*sqrt(N_M)
z_F<-(F_marginal+err_beta_F_scale)*sqrt(N_F)  

setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/moderate/",h_c[i]))
dir.create("causal")
save(select_cd,file = paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/moderate/",h_c[i],"/causal/data_",rep,".RData"))

#univariate sex-dimorphic analysis
dir.create("sd")
sd=(z_M/sqrt(N_M-1)-z_F/sqrt(N_F-1))/sqrt(1/(N_M-1)+1/(N_F-1))
log_p_value <- pnorm(abs(sd), lower.tail = FALSE, log.p = TRUE)
sdp <- 2 * exp(log_p_value)
snp=sub("_[^_]*$", "", rownames(sdp))
sdr=data.frame(snp,sdp)
save(sdr, file = paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/moderate/",h_c[i],"/sd/data_",rep,".RData"))

#stepwise regression
ref_M_bim<-fread(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/MALE/chr_",chr,".bim"))
data=ref_M_bim[match(snp,ref_M_bim$V2),][,c(2,5,6)]
colnames(data)=c("SNP","A1","A2")
N=N_M+N_F
data1=data
data1$b=(sqrt(N_M-1)*z_M+sqrt(N_F-1)*z_F)/(N-1)
data1$se=1/sqrt(N-1)
log_p_value <- pnorm(abs(data1$b/data1$se), lower.tail = FALSE, log.p = TRUE)
data1$p=2 * exp(log_p_value)
data1$N=N
data2=data
data2$b=(-sqrt(N_M-1)*z_M+sqrt(N_F-1)*z_F)/(N-1)
data2$se=1/sqrt(N-1)
log_p_value <- pnorm(abs(data2$b/data2$se), lower.tail = FALSE, log.p = TRUE)
data2$p=2 * exp(log_p_value)
data2$N=N
data2$SNP=paste0("S*",data2$SNP)
data=rbind(data1,data2)
N=N_M+N_F
V1=(N_M-1)*LD_M+(N_F-1)*LD_F
V2=(-1)*(N_M-1)*LD_M+(N_F-1)*LD_F
LD=rbind(cbind(V1,V2),cbind(V2,V1))/(N-1)
res=sdstepwise(data, LD, p_cutoff=5e-8, collinear=0.9,ncore=10)
dir.create("stepwise")
write.table(res,paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/moderate/",h_c[i],"/stepwise/res_",rep,".txt"), quote=FALSE, row.names=FALSE)

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

dir.create("susie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/moderate/",h_c[i],"/susie"))
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

#SuSiE-modify
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/moderate/",h_c[i]))
dir.create("susieadd")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/moderate/",h_c[i],"/susieadd"))
dir.create("susieresult")
dir.create("susiePIP")

sd=(z_M/sqrt(N_M-1)-z_F/sqrt(N_F-1))/sqrt(1/(N_M-1)+1/(N_F-1))
snp=sub("_[^_]*$", "", rownames(sd))
N=N_M+N_F
library(susieR)
R=((N_M-1)*LD_M+(N_F-1)*LD_F)/(N-1)
fit = susie_rss(c(sd),R,N)
save(fit,file = paste0("susieresult/region_",rep,".RData"))
PIP=fit$pip
data=data.frame(snp,PIP)
write.table(data,paste0("susiePIP/region_",rep,".txt"),quote=FALSE)

#MESuSiE-modify
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/moderate/",h_c[i]))
dir.create("mesusie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/moderate/",h_c[i],"/mesusie"))
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
sdp <- 2 * exp(log_p_value)
sdr=data.frame(snp,sdp)
save(sdr, file = paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/high/",h_c[i],"/mesusie/data_",rep,".RData"))

#SharePro_gxe
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/moderate/",h_c[i]))
dir.create("sharepro05")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/moderate/",h_c[i],"/sharepro05"))
SNP=sub("_[^_]*$", "", rownames(z_M))
N=N_M
Z=z_M
BETA=Z/sqrt(N_M-1)
SE=1/sqrt(N_M-1)
data1=data.frame(SNP,N,BETA,SE)

Z=z_F
BETA=Z/sqrt(N_F-1)
SE=1/sqrt(N_F-1)
N=N_F
data2=data.frame(SNP,N,BETA,SE)

write.table(data1,paste0("data1_",rep,".txt"),quote=F,row.names=F)
write.table(data2,paste0("data2_",rep,".txt"),quote=F,row.names=F)

write.table((LD_M+LD_F)/2,paste0("LD_",rep,".ld"),quote=F,row.names=F,col.names=F)

system(paste0("python /net/mulan/home/luliuu/SharePro_gxe/src/SharePro/sharepro_gxe.py --z data1_",rep,".txt data2_",rep,".txt --pthres 0.5 --ld LD_",rep,".ld --save res_",rep))
system(paste0("rm data1_",rep,".txt"))
system(paste0("rm data2_",rep,".txt"))
system(paste0("rm LD_",rep,".ld"))
}

#low LD intensity
set.seed(1217-rep)
select_cd=which.min(colSums(LD_M^2))
setwd("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore")
dir.create("low")

#set the parameter
h_c=c(1e-4,5e-4,1e-3)
h_d=5e-4
n_causal=1
ref_N_M=50000
ref_N_F=50000
N_M=150000
N_F=150000
rho=0.5

for(i in 1:3){
setwd("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low")
dir.create(as.character(h_c[i]))
h_d=h_c[i]
beta_c=rep(0,dim(LD_M)[1])
beta_d=rep(0,dim(LD_M)[1])

beta_shared<-rmvnorm(length(select_cd),mean=rep(0,2),sigma=matrix(c(h_c[i],rho*sqrt(h_d*h_c[i]),rho*sqrt(h_d*h_c[i]),h_d)/n_causal,ncol=2,nrow=2))
beta_c[select_cd]=beta_shared[,1]
beta_d[select_cd]<-beta_shared[,2]

M_marginal<-LD_M%*%(beta_c-beta_d)
F_marginal<-LD_F%*%(beta_c+beta_d)

y_null_M<-rnorm(ref_N_M,0,sqrt(1-h_c[i]-h_d))
err_beta_M<-t(ref_M_spe)%*%y_null_M/ref_N_M
y_null_F<-rnorm(ref_N_F,0,sqrt(1-h_c[i]-h_d))
err_beta_F<-t(ref_F_spe)%*%y_null_F/ref_N_F

err_beta_M_scale<-sqrt(ref_N_M/N_M)*err_beta_M
err_beta_F_scale<-sqrt(ref_N_F/N_F)*err_beta_F
      
z_M<-(M_marginal+err_beta_M_scale)*sqrt(N_M)
z_F<-(F_marginal+err_beta_F_scale)*sqrt(N_F)  

setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low/",h_c[i]))
dir.create("causal")
save(select_cd,file = paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low/",h_c[i],"/causal/data_",rep,".RData"))

#univariate sex-dimorphic analysis
dir.create("sd")
sd=(z_M/sqrt(N_M-1)-z_F/sqrt(N_F-1))/sqrt(1/(N_M-1)+1/(N_F-1))
log_p_value <- pnorm(abs(sd), lower.tail = FALSE, log.p = TRUE)
sdp <- 2 * exp(log_p_value)
snp=sub("_[^_]*$", "", rownames(sdp))
sdr=data.frame(snp,sdp)
save(sdr, file = paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low/",h_c[i],"/sd/data_",rep,".RData"))

#stepwise regression
ref_M_bim<-fread(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex/ref/MALE/chr_",chr,".bim"))
data=ref_M_bim[match(snp,ref_M_bim$V2),][,c(2,5,6)]
colnames(data)=c("SNP","A1","A2")
N=N_M+N_F
data1=data
data1$b=(sqrt(N_M-1)*z_M+sqrt(N_F-1)*z_F)/(N-1)
data1$se=1/sqrt(N-1)
log_p_value <- pnorm(abs(data1$b/data1$se), lower.tail = FALSE, log.p = TRUE)
data1$p=2 * exp(log_p_value)
data1$N=N
data2=data
data2$b=(-sqrt(N_M-1)*z_M+sqrt(N_F-1)*z_F)/(N-1)
data2$se=1/sqrt(N-1)
log_p_value <- pnorm(abs(data2$b/data2$se), lower.tail = FALSE, log.p = TRUE)
data2$p=2 * exp(log_p_value)
data2$N=N
data2$SNP=paste0("S*",data2$SNP)
data=rbind(data1,data2)
N=N_M+N_F
V1=(N_M-1)*LD_M+(N_F-1)*LD_F
V2=(-1)*(N_M-1)*LD_M+(N_F-1)*LD_F
LD=rbind(cbind(V1,V2),cbind(V2,V1))/(N-1)
res=sdstepwise(data, LD, p_cutoff=5e-8, collinear=0.9,ncore=10)
dir.create("stepwise")
write.table(res,paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low/",h_c[i],"/stepwise/res_",rep,".txt"), quote=FALSE, row.names=FALSE)

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

dir.create("susie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low/",h_c[i],"/susie"))
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

#SuSiE-modify
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low/",h_c[i]))
dir.create("susieadd")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low/",h_c[i],"/susieadd"))
dir.create("susieresult")
dir.create("susiePIP")

sd=(z_M/sqrt(N_M-1)-z_F/sqrt(N_F-1))/sqrt(1/(N_M-1)+1/(N_F-1))
snp=sub("_[^_]*$", "", rownames(sd))
N=N_M+N_F
library(susieR)
R=((N_M-1)*LD_M+(N_F-1)*LD_F)/(N-1)
fit = susie_rss(c(sd),R,N)
save(fit,file = paste0("susieresult/region_",rep,".RData"))
PIP=fit$pip
data=data.frame(snp,PIP)
write.table(data,paste0("susiePIP/region_",rep,".txt"),quote=FALSE)

#MESuSiE-modify
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low/",h_c[i]))
dir.create("mesusie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low/",h_c[i],"/mesusie"))
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
sdp <- 2 * exp(log_p_value)
sdr=data.frame(snp,sdp)
save(sdr, file = paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low/",h_c[i],"/mesusie/data_",rep,".RData"))

#SharePro_gxe
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low/",h_c[i]))
dir.create("sharepro05")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/LDscore/low/",h_c[i],"/sharepro05"))
SNP=sub("_[^_]*$", "", rownames(z_M))
N=N_M
Z=z_M
BETA=Z/sqrt(N_M-1)
SE=1/sqrt(N_M-1)
data1=data.frame(SNP,N,BETA,SE)

Z=z_F
BETA=Z/sqrt(N_F-1)
SE=1/sqrt(N_F-1)
N=N_F
data2=data.frame(SNP,N,BETA,SE)

write.table(data1,paste0("data1_",rep,".txt"),quote=F,row.names=F)
write.table(data2,paste0("data2_",rep,".txt"),quote=F,row.names=F)

write.table((LD_M+LD_F)/2,paste0("LD_",rep,".ld"),quote=F,row.names=F,col.names=F)

system(paste0("python /net/mulan/home/luliuu/SharePro_gxe/src/SharePro/sharepro_gxe.py --z data1_",rep,".txt data2_",rep,".txt --pthres 0.5 --ld LD_",rep,".ld --save res_",rep))
system(paste0("rm data1_",rep,".txt"))
system(paste0("rm data2_",rep,".txt"))
system(paste0("rm LD_",rep,".ld"))
}
