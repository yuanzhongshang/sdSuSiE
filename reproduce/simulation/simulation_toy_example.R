args <- commandArgs(TRUE)
rep <- as.integer(args[1])

set.seed(1217-rep)
library(data.table)
library(sdSuSiE)
library(mvtnorm)

setwd("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/simple")

#set the parameter
h_c=c(1e-4,1e-5,1e-6,1e-7)
n_causal=1
ref_N_M=50000
ref_N_F=50000
N_M=150000
N_F=150000
rho=0.5
snp=c("X1","X2","X3")

LD=matrix(c(1,1,0,1,1,0,0,0,1),3,3)
ref_M_spe=rmvnorm(50000,mean=rep(0,3),sigma=LD)
ref_F_spe=rmvnorm(50000,mean=rep(0,3),sigma=LD)
ref_M_spe=scale(ref_M_spe)
ref_F_spe=scale(ref_F_spe)
LD_M=cor(ref_M_spe)
LD_F=cor(ref_F_spe)

for(i in 1:4){
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/simple/"))
dir.create(as.character(h_c[i]))
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/simple/",h_c[i]))

h_d=h_c[i]
beta_c=rep(0,dim(LD_M)[1])
beta_d=rep(0,dim(LD_M)[1])

beta_shared<-rmvnorm(1,mean=rep(0,2),sigma=matrix(c(h_c[i],rho*sqrt(h_d*h_c[i]),rho*sqrt(h_d*h_c[i]),h_d)/n_causal,ncol=2,nrow=2))
beta_c[1]=beta_shared[,1]
beta_d[1]<-beta_shared[,2]

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

dir.create("1")
dir.create("0")

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

setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/simple/",h_c[i],"/1"))
dir.create("susie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/simple/",h_c[i],"/1/susie"))
dir.create("susieresult")
dir.create("susieresultL_10")
dir.create("susiePIP")
summary_list=list()
summary_list[[1]]=data1[1:2,]
summary_list[[2]]=data2[1:2,]
names(summary_list)=c("Male","Female")

mat_list=list()
mat_list[[1]]=LD_M[1:2,1:2]
mat_list[[2]]=LD_F[1:2,1:2]
names(mat_list)=c("Male","Female")

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

setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/simple/",h_c[i],"/0"))
dir.create("susie")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/simple/",h_c[i],"/0/susie"))
dir.create("susieresult")
dir.create("susieresultL_10")
dir.create("susiePIP")
summary_list=list()
summary_list[[1]]=data1[c(1,3),]
summary_list[[2]]=data2[c(1,3),]
names(summary_list)=c("Male","Female")

mat_list=list()
mat_list[[1]]=LD_M[c(1,3),c(1,3)]
mat_list[[2]]=LD_F[c(1,3),c(1,3)]
names(mat_list)=c("Male","Female")

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
  
#susie for the SDE
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/simple/",h_c[i],"/1"))
dir.create("susieadd")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/simple/",h_c[i],"/1/susieadd"))
dir.create("susieresult")
dir.create("susiePIP")

sd=(z_M/sqrt(N_M-1)-z_F/sqrt(N_F-1))/sqrt(1/(N_M-1)+1/(N_F-1))
N=N_M+N_F
library(susieR)
R=((N_M-1)*LD_M+(N_F-1)*LD_F)/(N-1)
fit = susie_rss(c(sd)[1:2],R[1:2,1:2],N)
save(fit,file = paste0("susieresult/region_",rep,".RData"))
PIP=fit$pip
data=data.frame(snp[1:2],PIP)
write.table(data,paste0("susiePIP/region_",rep,".txt"),quote=FALSE)

setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/simple/",h_c[i],"/0"))
dir.create("susieadd")
setwd(paste0("/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/20250404/simple/",h_c[i],"/0/susieadd"))
dir.create("susieresult")
dir.create("susiePIP")
fit = susie_rss(c(sd)[c(1,3)],R[c(1,3),c(1,3)],N)
save(fit,file = paste0("susieresult/region_",rep,".RData"))
PIP=fit$pip
data=data.frame(snp[c(1,3)],PIP)
write.table(data,paste0("susiePIP/region_",rep,".txt"),quote=FALSE)

}
