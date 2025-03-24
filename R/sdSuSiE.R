#' @title The main function for sex-dimorphic mapping with sum of the single effects model (sdSuSiE).
#' @description sdSuSiE conduct sex-stratified GWAS fine-mapping to identify potentially sex-dimorphic causal SNPs associated with a complex trait of interest. 
#' @param R_mat_list A list of length 2, where:
#' - The first element is the correlation matrix for males.
#' - The second element is the correlation matrix for females.
#' The column names of each correlation matrix must match the order of SNP names in the corresponding element of `summary_stat_list`.
#'
#' @param summary_stat_list A list of length 2, where:
#' - The first element contains the summary statistics for males.
#' - The second element contains the summary statistics for females.
#' Each summary statistics element must include the following columns: SNP, Beta, Se, Z, and N. 
#' The order of SNPs in each summary statistics element must match the order of the corresponding correlation matrix in `R_mat_list`. 
#' The function `sdSuSiE` reconstructs sufficient statistics of each SNP for each sex using either:
#' - Marginal z-scores (Z) and sample sizes (N), or
#' - Marginal effect sizes (Beta) and standard errors (Se).
#'
#' @param L The maximum number of non-zero effects assumed within the region (default: 10).
#' 
#' @param residual_variance A numeric vector (default: `c(1, 1)`) specifying the residual variance for males and females.
#' 
#' @param prior_weights A numeric vector of length p specifying the prior probability that each SNP has a non-zero effect. 
#' By default, prior weights are assumed to be equal for all SNPs.
#'
#' @param config_weights A numeric vector of length 3 specifying the prior probabilities for different causal configurations:
#' - The first element represents the probability of a SNP having a common effect.
#' - The second element represents the probability of a SNP having a sex-dimorphic effect.
#' - The third element represents the probability of a SNP having both common and sex-dimorphic effects.
#' The three elements must sum to 1. The default prior values are `-1 + sqrt(2)`, `-1 + sqrt(2)`, and `(-1 + sqrt(2))^2`, respectively.
#' 
#' @param config_names A character vector of length 3 specifying the names for different causal configurations. The default prior values are `Common`, `Sex-dimorphic`, and `Both`, respectively.
#'
#' @param estimate_residual_variance A logical value indicating whether the residual variance should be estimated (`TRUE`) or fixed (`FALSE`, default).
#'   
#' @param cor_method A string specifying the method to compute the 95% credible set based on correlations among SNPs in the region. Options include:
#' - `"min_abs_corr"`: Minimum absolute correlation (default).
#' - `"mean_abs_corr"`: Mean absolute correlation.
#' - `"median_abs_corr"`: Median absolute correlation.
#' 
#' @param cor_threshold A numeric value specifying the correlation threshold for `cor_method` (default: 0.5).
#' 
#' @param max_iter An integer specifying the maximum number of iterations to perform (default: 100).
#' 
#' @return An R6 object with PIPs, credible sets, and other features of the sex-dimorphic fine-mapping result. 
#' \item{PIP}{A vector where each element represents the PIP of a SNP having a non-zero effect.}
#' 
#' \item{PIP_config}{A matrix where each column corresponds to the PIP of a causal configuration. The columns represent the PIP of a SNP having a common effect, a sex-dimorphic effect, or both.}
#' 
#' \item{alpha}{A list of length L, where each element is a matrix (p by the number of causal configurations) of PIPs. Each column represents the PIP for a specific scenario.}
#' 
#' \item{KL}{A vector containing the single-effect Kullback–Leibler (KL) divergences.}
#' 
#' \item{sigma2}{A numeric vector of estimated residual variance.}
#' 
#' \item{V}{A list of length L, where each element is a 2 by 2 matrix representing the prior variance estimate.}
#' 
#' \item{ELBO}{A vector storing the evidence lower bound achieved at each iteration of the model-fitting algorithm, which aims to maximize the ELBO.}
#' 
#' \item{Credible_Set}{The estimated credible sets for fine-mapping.}
#' 
#' \item{posterior_mean}{A matrix where the first column contains the posterior mean for \beta_c (common effects) and the second column contains the posterior mean for \eqn{\beta_d} (sex-dimorphic effects).}
#' 
#' \item{posterior_variance}{A matrix with the first column being the variance for \eqn{\beta_c}, the second column being the variance for \eqn{\beta_d}, and the third column being the covariance between \eqn{\beta_c} and \eqn{\beta_d}.}
#'
#' @import R6 Rcpp RcppArmadillo nloptr
#' @export
#' @examples
#' library(sdSuSiE)
#' data(GWAS_M)
#' data(GWAS_F)
#' data(LD_M)
#' data(LD_F)
#' summary_list = list()
#' summary_list[[1]] = GWAS_M
#' summary_list[[2]] = GWAS_F
#' names(summary_list) = c("Male","Female")
#' LDmat_list = list()
#' LDmat_list[[1]] = LD_M
#' LDmat_list[[2]] = LD_F
#' names(LDmat_list) = c("Male","Female")
#' res <- sdSuSiE(LDmat_list, summary_list)

sdSuSiE<-function(R_mat_list, summary_stat_list, L = 10, residual_variance = NULL, prior_weights = NULL, config_weights = NULL, config_names = NULL, estimate_residual_variance = FALSE, cor_method = "min.abs.corr", cor_threshold = 0.5, max_iter = 100){
  
  cat("*************************************************************\n
  Sex-dimorphic mapping with sum of the single effects model (sdSuSiE)          \n
   Visit https://xiangzhou.github.io/software/ For Update            \n
            (C) 2025 Lu Liu, Xiang Zhou                        \n
              GNU General Public License                          \n
*************************************************************") 
  
  time_start <- Sys.time()
  cat("\n# Start data processing for sufficient statistics \n")
  sdSuSiEData_obj <- sdSuSiEData$new(R_mat_list, summary_stat_list)
  
  n_snp = nrow(summary_stat_list[[1]])
  n_sex = length(summary_stat_list)
  
  if(is.null(prior_weights)){
	prior_weights = rep(1/n_snp, n_snp)
   }
   
  if(is.null(config_weights)){
    base_fac_ratio = (-1+sqrt(2))/(-1+sqrt(2))^2
	base_fac = 1/Reduce("+", lapply(seq(1,n_sex), function(x)choose(n_sex, x) * base_fac_ratio^(n_sex - x)))
	config_weights = unlist(lapply(seq(1,n_sex), function(x)rep(base_fac * base_fac_ratio^(n_sex - x), choose(n_sex, x))))
   }
   used_weights = kronecker(prior_weights, t(config_weights), FUN = "*")

  if(is.null(residual_variance)){
    residual_variance = rep(1, n_sex)
  }
  
  if(is.null(config_names)){
	config_names = c("Common", "Sex-dimorphic", "Both")
  }
  
  cat("# Create sdSuSiE object \n")
  sdSuSiEObject_obj <- sdSuSiEObject$new(n_snp, L, n_sex, residual_variance, used_weights, estimate_residual_variance, max_iter, config_names)
  cat("# Start data analysis \n")
  
  n_iter = 0
  for (iter in 1:max_iter) {
    
    comp_residual <- sdSuSiEObject_obj$compute_residual(sdSuSiEData_obj, sdSuSiEObject_obj)
    
    
    for (l_index in seq(1,L,1)) {
      
      comp_residual = comp_residual + sdSuSiEObject_obj$Xr[[l_index]]
      
      SER_res <- single_effect_regression(comp_residual, sdSuSiEData_obj$XtX.diag, sdSuSiEObject_obj, l_index) 
      
      sdSuSiEObject_obj$par_update(SER_res, l_index)
      	  
      sdSuSiEObject_obj$compute_KL(SER_res, sdSuSiEObject_obj$compute_SER_posterior_loglik(sdSuSiEData_obj, comp_residual, SER_res$b1b2), l_index)
      
      sdSuSiEObject_obj$compute_Xr(sdSuSiEData_obj, SER_res$b1b2$EB1, l_index)
      
      comp_residual = comp_residual - sdSuSiEObject_obj$Xr[[l_index]]
         
    }

    updated_sigma2 = sdSuSiEObject_obj$update_residual_variance(sdSuSiEData_obj, iter)
    
    if((sdSuSiEObject_obj$ELBO[iter+1] - sdSuSiEObject_obj$ELBO[iter]) < 0.001){
      break
    }
	
    if(sdSuSiEObject_obj$estimate_residual_variance == TRUE){
      sdSuSiEObject_obj$sigma2 =  updated_sigma2
    }
	
    n_iter = n_iter + 1
    #check convergence and update sigma2
  }

  cat("\n# Data analysis is done, and now generates result \n\n")
   
  sdSuSiEObject_obj$get_result(sdSuSiE_get_cs(sdSuSiEObject_obj, R_mat_list, cor_method = cor_method, cor_threshold = cor_threshold), sdSuSiE_get_pip_either(sdSuSiEObject_obj), sdSuSiE_get_pip_config(sdSuSiEObject_obj))
  sdSuSiEObject_obj$sdSuSiE_summary(sdSuSiEData_obj)

  time_end <- Sys.time()
  cat(c("\n# Total time used for the analysis:", paste0(round(as.numeric(difftime(time_end, time_start, units = c("mins"))), 2), " mins\n")))
  return(sdSuSiEObject_obj)
}


#' sdSuSiE Object
#'
#' This object is designed to manage and perform sdSuSiE analysis.
#'
#' @name sdSuSiEObject
#' @description An R6 object for handling the sdSuSiE analysis and storing its results.
#' 
#' @field name_config A character vector containing the names of causal configurations.
#' @field column_config A list specifying the possible configurations for columns.
#' @field alpha A list of alpha matrices, where each matrix corresponds to posterior inclusion probabilities for each component.
#' @field mu1 A list of mu1 lists, where each list contains the posterior means for each causal configuration for each component.
#' @field mu2 A list of mu2 lists, where each list contains the posterior variances for each causal configuration for each component.
#' @field EB1 A list of EB1 matrices, where each matrix represents the posterior means of \eqn{\beta_c} and \eqn{\beta_d} for each component.
#' @field EB2 A list of EB2 matrices, where each matrix represents the posterior variances for each component. The first column is the variance for \eqn{\beta_c}, the second column is the variance for \eqn{\beta_d}, and the third column is the covariance between \eqn{\beta_c} and \eqn{\beta_d}.
#' @field Xr A list of Xr matrices, representing the genetic part for each component.
#' @field KL A numeric vector containing the Kullback-Leibler (KL) divergence values for each component.
#' @field lbf A list of lbf matrices, where each matrix represents the log-Bayes factors for causal configuration for each component.
#' @field ELBO A numeric vector of Evidence Lower Bound (ELBO) values for each iteration, used to assess model convergence.
#' @field sigma2 A numeric vector of residual variances.
#' @field V A list of V matrices, representing the prior variance for each component.
#' @field pi A numeric vector of prior weights for each component.
#' @field estimate_residual_variance A logical value indicating whether the residual variance is estimated (\code{TRUE}) or fixed (\code{FALSE}).
#' @field L An integer specifying the number of causal components in the model.
#' @field nSNP An integer representing the total number of SNPs analyzed.
#' @field nsex An integer representing the number of sexes analyzed (typically 2: male and female).
#' @field Credible_Set A list of credible sets, each representing SNP sets with high posterior probabilities of causality.
#' @field PIP A numeric vector of Posterior Inclusion Probabilities (PIP) for each SNP.
#' @field PIP_config A matrix of PIP values for each SNP across different configurations.
#' @field posterior_mean A matrix containing the posterior means for \eqn{\beta_c} (common effects) and \eqn{\beta_d} (sex-dimorphic effects).
#' @field posterior_variance A matrix with the first column being the variance for \eqn{\beta_c}, the second column being the variance for \eqn{\beta_d}, and the third column being the covariance between \eqn{\beta_c} and \eqn{\beta_d}.
#'
#' @import R6
#' @export

sdSuSiEObject <- R6::R6Class("sdSuSiEObject", public = list(
  #' @description Initialize the sdSuSiE object.
  #' @param p Integer, number of SNPs.
  #' @param L Integer, number of causal components.
  #' @param N_sex Integer, number of sexes analyzed (typically 2: male and female).
  #' @param residual_variance Numeric vector of initial residual variances.
  #' @param prior_weights Numeric vector of prior weights.
  #' @param estimate_residual_variance Logical, if TRUE, estimate the residual variance.
  #' @param max_iter Integer, maximum number of iterations.
  #' @param name_vec Character vector, names for configurations.
  #' @return A new `sdSuSieObject` object.
  initialize = function(p, L, N_sex, residual_variance, prior_weights, estimate_residual_variance, max_iter, name_vec){
    
	self$name_config = name_vec
	self$column_config = Reduce(append, lapply(seq(1, N_sex), function(x){
		poss_config = combn(1:N_sex, x)
		lapply(1:ncol(poss_config), function(y)return(matrix(poss_config[,y])))
	}))
	
	self$alpha = rep(list(matrix(0, nrow = p, ncol = N_sex)), L)
    self$mu1 = rep(list(lapply(self$column_config, function(x){matrix(0, nrow = p,ncol = length(x))})), L)
    self$mu2 = rep(list(lapply(c(1,1,3), function(x){matrix(0, nrow = p, ncol = x)})), L)
	self$EB1 = rep(list(matrix(0, nrow = p, ncol = N_sex)), L)
	self$EB2 = rep(list(matrix(0, nrow = p, ncol = 3)), L)
	
    self$Xr = rep(list(matrix(0,nrow = p, ncol = N_sex)), L)
    self$KL = rep(as.numeric(NA), L)
    self$lbf = vector("list", L)
    self$ELBO = rep(NA, max_iter)
    self$ELBO[1] = -Inf
    self$sigma2 = residual_variance
    self$V = rep(list(matrix(0, ncol = N_sex, nrow = N_sex)), L)
    self$pi = prior_weights
    self$estimate_residual_variance = estimate_residual_variance   
    self$L = L
    self$nSNP = p
    self$nsex = N_sex
    
    self$Credible_Set = list()
    self$PIP = rep(as.numeric(NA), p)
	self$PIP_config = matrix(as.numeric(NA), p, 2^N_sex-1)
	self$posterior_mean = matrix(as.numeric(NA), p, 2)
	self$posterior_variance = matrix(as.numeric(NA), p, 3)
  },    

  #' @param sdSuSiE_Data List, sdSuSiE data object.
  #' @param sdSuSiE_Obj, current sdSuSiE object.
  #' @return Matrix of computed residuals.
  compute_residual = function (sdSuSiE_Data, sdSuSiE_Obj){
    residual <- Reduce(cbind, sdSuSiE_Data$Xty_list) - Reduce("+", sdSuSiE_Obj$Xr)
    return(residual)
  },

  #' @param sdSuSiE_Data List, sdSuSiE data object.
  #' @param b1b2 Matrix, current estimate of b.
  #' @param l_index Integer, index of the current component.
  #' @return Updated Xr matrix for the specified component.
  compute_Xr = function(sdSuSiE_Data, b1b2, l_index){
	c1 = sdSuSiE_Data$XtX_list[[1]] %*% (b1b2[,1] - b1b2[,2])
	c2 = sdSuSiE_Data$XtX_list[[2]] %*% (b1b2[,1] + b1b2[,2])
    self$Xr[[l_index]] = cbind(c1,c2)
  },
  
  #' @param SER_res List, results from the SER analysis.
  #' @param l_index Integer, index for the current component.
  #' @return Updated `sdSuSiEObject` object containing the analysis results and model parameters.
  par_update = function(SER_res, l_index){
    self$alpha[[l_index]] <- SER_res$alpha
    self$mu1[[l_index]] <- SER_res$mu1_multi
    self$mu2[[l_index]] <- SER_res$mu2_multi
    self$lbf[[l_index]] <- SER_res$lbf_multi
    self$V[[l_index]] <- SER_res$V
	self$EB1[[l_index]] <- SER_res$b1b2$EB1
	self$EB2[[l_index]] <- SER_res$b1b2$EB2
  },

  #' @param sdSuSiE_Data List, sdSuSiE data object.
  #' @param comp_residual Matrix, residual values for the current component.
  #' @param b1b2 List, containing the EB1 and EB2 matrices.
  #' @return Log-likelihood value for the SER posterior.
   compute_SER_posterior_loglik = function(sdSuSiE_Data, comp_residual, b1b2){
    return(
      sum(-0.5/self$sigma2[1]*(-2*comp_residual[,1]*(b1b2$EB1[,1]-b1b2$EB1[,2])+sdSuSiE_Data$XtX.diag[[1]]*(b1b2$EB2[,1]+b1b2$EB2[,2]-2*b1b2$EB2[,3])))+
	  sum(-0.5/self$sigma2[2]*(-2*comp_residual[,2]*(b1b2$EB1[,1]+b1b2$EB1[,2])+sdSuSiE_Data$XtX.diag[[2]]*(b1b2$EB2[,1]+b1b2$EB2[,2]+2*b1b2$EB2[,3])))
    )
  },
  
  #' @param SER_res List, results from the SER analysis.
  #' @param value Numeric, value to compute KL divergence.
  #' @param l_index Integer, index of the current component.
  #' @return Updated KL divergence value for the specified component.  
  compute_KL = function(SER_res, value, l_index){  
    self$KL[l_index] = -SER_res$loglik + value    
  },

  #' @param sdSuSiE_Data List, sdSuSiE data object.
  #' @param niter Integer, the current iteration number.
  #' @return A numeric vector with the updated residual variances.  
  update_residual_variance = function(sdSuSiE_Data, niter){    
    B_1 = lapply(1:self$nsex, function(y)Reduce(cbind, lapply(self$EB1, function(x)x[,y])))   
    BXXB = list()
	BXXB[[1]] = sum((t(B_1[[1]]-B_1[[2]])%*%sdSuSiE_Data$XtX_list[[1]])*t(B_1[[1]]-B_1[[2]]))
	BXXB[[2]] = sum((t(B_1[[1]]+B_1[[2]])%*%sdSuSiE_Data$XtX_list[[2]])*t(B_1[[1]]+B_1[[2]]))
			
    betabar = Reduce("+", self$EB1)
    
	BbarXXBbar = list()
	BbarXXBbar[[1]] = sum((betabar[,1]-betabar[,2])*(sdSuSiE_Data$XtX_list[[1]]%*%(betabar[,1]-betabar[,2])))
	BbarXXBbar[[2]] = sum((betabar[,1]+betabar[,2])*(sdSuSiE_Data$XtX_list[[2]]%*%(betabar[,1]+betabar[,2])))

	BbarXty = list()
	BbarXty[[1]] = 2*sum((betabar[,1]-betabar[,2])*sdSuSiE_Data$Xty_list[[1]])
	BbarXty[[2]] = 2*sum((betabar[,1]+betabar[,2])*sdSuSiE_Data$Xty_list[[2]])
	
    B_2 = lapply(1:3, function(y)Reduce(cbind, lapply(self$EB2, function(x)x[,y])))
	XXB = list()
	XXB[[1]] = sum(sdSuSiE_Data$XtX.diag[[1]]*(B_2[[1]]+B_2[[2]]-2*B_2[[3]]))
	XXB[[2]] = sum(sdSuSiE_Data$XtX.diag[[2]]*(B_2[[1]]+B_2[[2]]+2*B_2[[3]]))

	updated_sigma = lapply(1:self$nsex, function(x){
      ((sdSuSiE_Data$yty_list[[x]]-BbarXty[[x]]+BbarXXBbar[[x]])+(XXB[[x]]-BXXB[[x]]))/sdSuSiE_Data$N_list[[x]]
    })
	
    self$ELBO[niter+1] = Reduce(sum, lapply(1:self$nsex, function(x){
      -sdSuSiE_Data$N_list[[x]]/2*log(2*pi*self$sigma2[x])-(1/(2*self$sigma2[x]))*updated_sigma[[x]]*sdSuSiE_Data$N_list[[x]]
    })) - sum(self$KL)
    
    return(unlist(updated_sigma))
  },

  #' @param cs List, credible sets for effects.
  #' @param pip Numeric vector, posterior inclusion probabilities for each SNP.
  #' @param pip_config Matrix, PIP values for each SNP across  
  get_result = function(cs, pip, pip_config){
    self$Credible_Set = cs
    self$PIP = pip
	self$PIP_config = pip_config
	colnames(self$PIP_config) = self$name_config
	self$posterior_mean <- Reduce("+", self$EB1)
	colnames(self$posterior_mean) = c("beta_c", "beta_d")
	var_c = 0
	var_d = 0
	var_cd = 0
	for(j in 1:length(self$alpha)){
	  var_c=var_c+self$alpha[[j]][,1]*(self$mu2[[j]][[1]]-self$mu1[[j]][[1]]^2)+self$alpha[[j]][,3]*(self$mu2[[j]][[3]][,1]-self$mu1[[j]][[3]][,1]^2)
	  var_d=var_d+self$alpha[[j]][,2]*(self$mu2[[j]][[2]]-self$mu1[[j]][[2]]^2)+self$alpha[[j]][,3]*(self$mu2[[j]][[3]][,2]-self$mu1[[j]][[3]][,2]^2)
	  var_cd=var_cd+self$alpha[[j]][,3]*(self$mu2[[j]][[3]][,3]-self$mu1[[j]][[3]][,2]*self$mu1[[j]][[3]][,1])
	}
	self$posterior_variance = data.frame(var_c, var_d, var_cd)
  },

  #' @param sdSuSiE_Data List, sdSuSiE data object.
  sdSuSiE_summary = function(sdSuSiE_Data){
    cat(c(paste0("Potential causal SNPs with common effect (PIPc > 0.5): "), sdSuSiE_Data$Summary_Stat[[1]]$SNP[which(self$PIP_config[,1]>0.5)], "\n\n"))
    cat(c(paste0("Potential causal SNPs with sex-dimorphic effect (PIPd > 0.5): "), sdSuSiE_Data$Summary_Stat[[1]]$SNP[which(self$PIP_config[,2]>0.5)], "\n\n"))
	cat(c(paste0("Potential causal SNPs with both common and sex-dimorphic effects (PIPb > 0.5): "), sdSuSiE_Data$Summary_Stat[[1]]$SNP[which(self$PIP_config[,3]>0.5)], "\n\n"))
	cat("Credible sets for effects: \n")
    print(self$Credible_Set)
    cat("\n Use sdSuSiE_Plot() for visualization")
  }
  
),lock_objects = F
)


#' sdSuSiEData
#' 
#' This R6 class is designed for managing and processing data relevant to the sdSuSiE approach.
#' The class includes methods to initialize and process XtX and XtY data structures.
#' 
#' @name sdSuSiEData
#' @description An R6 class for managing and processing data relevant to the sdSuSiE approach.
#'
#' @field R A list of correlation matrix. Each element represents a correlation matrix for a specific group (e.g., males and females).
#' @field Summary_Stat A list of summary statistics. Each element contains the required columns, such as SNP, Beta, Se, Z, and N, for a specific group.
#' @field var_y Variance of y. Default is 1.
#' @field Name_list A list of SNP names extracted from the correlation matrices in `R`.
#' @field N_sex Number of sexs in `X`.
#' @field XtX.diag A vector containing the diagonal elements of the XtX matrix.
#' @field XtX_list A list of processed XtX matrices for each group.
#' @field Xty_list A list of processed XtY vectors for each group.
#' @field N_list A list of median sample sizes (`N`) for each group, computed from `Summary_Stat`.
#' @field yty_list A list of precomputed yty values for each group.
#' 
#' @return An updated `sdSuSiEData` class is used for handling and processing data in the sdSuSiE method.
#' @import R6
#' @export

sdSuSiEData <- R6::R6Class("sdSuSiEData", public = list(
  #' @description Initialize the sdSuSiEData class
  #' @param X A list or data frame containing the input correlation matrix.
  #' @param Y A list or data frame containing the summary statistics.
  #' @param var_y A numeric value representing the variance of Y, default is 1.
  initialize = function(X, Y, var_y = 1){
    self$R <- X
    self$Summary_Stat <- Y
    self$var_y = var_y    
    self$Name_list <- as.list(names(self$R))
    names(self$Name_list) <- names(self$R)
    self$N_sex <- length(X)
    self$N_list <- lapply(self$Summary_Stat, function(x)median(x$N))
    self$yty_list <- lapply(self$N_list, function(x)return(self$var_y*(x-1)))
    
    self$XtX.diag <- self$XtX_diag(self$Summary_Stat, self$Name_list)
    self$XtX_list <- self$XtX_pro(self$R, self$XtX.diag, self$Name_list)
    self$Xty_list <- self$Xty_pro(self$Summary_Stat, self$XtX.diag, self$Name_list) ##diag(xtx)^*betahat
    return(self)
  },
	
  #' @description Compute XtX.diag for each summary statistic
  #' @param Summary_Stat A list containing the summary statistics.
  #' @param Name_list A list of names corresponding to the summary statistics.
  #' @return A list of computed XtX.diag values for each summary statistic.
  XtX_diag = function(Summary_Stat, Name_list){
    return( lapply(Name_list,function(x){
      R2 = (Summary_Stat[[x]]$Z^2)/(Summary_Stat[[x]]$Z^2+Summary_Stat[[x]]$N-2)
      sigma2 = self$var_y*(1-R2)*(Summary_Stat[[x]]$N-1)/(Summary_Stat[[x]]$N-2)
      return(sigma2/(Summary_Stat[[x]]$Se)^2)
    }))
  },
  
  #' @description Process XtX using the given R matrix and XtX.diag values
  #' @param R A list or data frame representing the correlation matrix.
  #' @param XtX.diag A list of XtX.diag values computed earlier.
  #' @param Name_list A list of names corresponding to the summary statistics.
  #' @return A list of processed XtX matrices for each summary statistic.
  XtX_pro = function(R, XtX.diag, Name_list){
    return(lapply(Name_list,function(x){
      return(diag(sqrt(XtX.diag[[x]]))%*%R[[x]]%*%diag(sqrt(XtX.diag[[x]])))
    }))
  },
  
  #' @description Process XtY (product of XtX.diag and Beta for each summary statistic)
  #' @param Summary_Stat A list containing the summary statistics.
  #' @param XtX.diag A list of XtX.diag values computed earlier.
  #' @param Name_list A list of names corresponding to the summary statistics.
  #' @return A list of computed XtY values for each summary statistic.
  Xty_pro = function(Summary_Stat, XtX.diag, Name_list){
    return(lapply(Name_list,function(x){
      XtX.diag[[x]]*Summary_Stat[[x]]$Beta
    }))
  }
),
lock_objects = F)

#' @title The function for single effect regression in sdSuSiE.
#' @description single_effect_regression conducts single effect regression for the specific component. 
#'
#' @param XtR A list of Xr matrices, representing the genetic part for the specific component.
#' @param XtX.diag A numeric vector containing the diagonal elements of the XtX matrix.
#' @param sdSuSiEObject_obj An R6 object for handling the sdSuSiE analysis and storing its results.
#' @param l_index A numeric value specifying the index of the component for which single-effect regression is performed.
#' 
#' @return A list with alphas, posteriors, likelihood, prior variance estimates, for single effect regression result for the specific component. 
#' \item{alpha}{A matrix of dimension \eqn{p \times k}, where \eqn{p} is the number of SNPs, and \eqn{k} is the number of causal configurations. Each column contains the PIPs for a specific causal scenario.}
#' \item{mu1_multi}{A list of posterior means for each causal configuration.}
#' \item{mu2_multi}{A list of posterior variances for each causal configuration.}
#' \item{lbf_multi}{A matrix containing log-Bayes factors for each causal configuration.}
#' \item{V}{A 2 by 2 matrix representing the prior variance estimates.}
#' \item{loglik}{The logarithm of the likelihood.} 
#' \item{b1b2}{A list containing EB1 and EB2 matrices.}
#' 
#' @import R6 nloptr Rcpp RcppArmadillo
#' @export

single_effect_regression <- function(XtR, XtX.diag, sdSuSiEObject_obj, l_index){
  column_config = sdSuSiEObject_obj$column_config
  N_sex = sdSuSiEObject_obj$nsex
  
  Xty_standardized = Reduce(cbind, lapply(1:N_sex, function(x){
    XtR[,x]/sdSuSiEObject_obj$sigma2[x]
  }))
  
  shat2 = Reduce(cbind, lapply(1:N_sex, function(x){
    sdSuSiEObject_obj$sigma2[x]/XtX.diag[[x]]
  }))
  
  betahat = shat2 * Xty_standardized

  opt_par<-pre_optim(N_sex, -30, 10)		
  if(N_sex == 2){		 
	update_V <- optim(opt_par$inital_par, fn = loglik_cpp, gr = NULL, betahat = betahat, shat2 = shat2, prior_weight = sdSuSiEObject_obj$pi, nsex = opt_par$nsex, diag_index = opt_par$diag_index, config_list = column_config, method = "L-BFGS-B", lower = opt_par$lower_bound, upper = opt_par$upper_bound)
	V_mat = vec_to_cov(update_V$par, opt_par$diag_index, opt_par$nsex)	 
  }else{
	intermediate_V <- nloptr(opt_par$inital_par, eval_f = loglik_cpp, eval_grad_f = NULL, lb = opt_par$lower_bound, ub = opt_par$upper_bound, betahat = betahat, shat2 = shat2, prior_weight = sdSuSiEObject_obj$pi, nsex = opt_par$nsex, diag_index = opt_par$diag_index, config_list = column_config, opts = list("algorithm" = "NLOPT_GN_DIRECT_L", "xtol_rel" = 1.0e-10))
	update_V <- nloptr(intermediate_V$solution, eval_f = loglik_cpp, eval_grad_f = NULL, lb = opt_par$lower_bound, ub = opt_par$upper_bound, betahat = betahat, shat2 = shat2, prior_weight = sdSuSiEObject_obj$pi, nsex = opt_par$nsex, diag_index = opt_par$diag_index, config_list = column_config, opts = list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1.0e-10))
	V_mat = vec_to_cov(update_V$solution, opt_par$diag_index, opt_par$nsex) 
  }

  reg_out <- lapply(column_config, function(x){
	if(length(x) == 1){
	  uni_reg(betahat, shat2, V_mat, x)
	}else if(length(x) > 1){
	  mvlmm_reg(betahat[,x], shat2[,x], V_mat[x,x]) 
	}
  })

  lbf <- Reduce(cbind, lapply(reg_out, function(x)x$lbf))
  lbf[is.na(lbf)] <- 0
  softmax_out <- compute_softmax(lbf, sdSuSiEObject_obj$pi)
  mu1_multi = lapply(reg_out, function(x)x$post_mean)
  mu2_multi = lapply(reg_out, function(x)x$post_mean2)
  b1b2 = compute_b1b2(softmax_out$alpha_wmulti, mu1_multi, mu2_multi, column_config, ncol(betahat), nrow(betahat))
	
  return(list(alpha = softmax_out$alpha_wmulti, mu1_multi = mu1_multi, mu2_multi = mu2_multi, lbf_multi = lbf, V = V_mat, loglik = softmax_out$loglik, b1b2 = b1b2))
}


#' @title The function prepares the input for the parameter space of the prior covariance matrix.
#'
#' @param nsex Number of sexes.
#' @param min_var Lower bound of the prior variance.
#' @param max_var Upper bound of the prior variance.
#' 
#' @return A list containing the initial parameters, lower bounds, upper bounds, diagonal element indices, and the number of sexes.
#' @export

pre_optim <- function(nsex, min_var, max_var) {
  
  # Compute the total number of elements in the covariance matrix
  total_ele_num <- (nsex + 1) * nsex / 2
  
  # Compute the indices of diagonal elements in the covariance matrix
  diag_ele_index <- do.call(rbind, lapply(seq(1, nsex), function(x) x * (x + 1) / 2))
  
  # Set initial parameters to zeros
  inital_par <- rep(0, total_ele_num)
  
  # Set default bounds for off-diagonal elements
  lower_bound <- rep(-0.9999, total_ele_num)
  upper_bound <- rep(0.9999, total_ele_num)
  
  # Override bounds for diagonal elements using provided min_var and max_var
  lower_bound[diag_ele_index] <- min_var
  upper_bound[diag_ele_index] <- max_var
  
  # Return a list containing the constructed parameters and bounds
  return(list(
    inital_par = inital_par,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    diag_index = diag_ele_index,
    nsex = nsex
  ))
}


#' @title The function transforms the estimated parameter vector into a covariance matrix.
#'
#' @param v_vec The estimated parameter vector.
#' @param diag_index Index of the variance in the parameter vector.
#' @param nsex Upper bound of the prior variance.
#' 
#' @return A constructed covariance matrix.
#' @export

vec_to_cov <- function(v_vec, diag_index, nsex) {
  
  # Initialize an empty correlation matrix filled with NAs
  cor_mat <- matrix(NA, ncol = nsex, nrow = nsex)
  
  # Fill the upper triangle (including diagonal) of the correlation matrix using the parameter vector
  cor_mat[upper.tri(cor_mat, diag = TRUE)] <- v_vec
  
  # Mirror the upper triangle values to the lower triangle to make the matrix symmetric
  cor_mat[lower.tri(cor_mat)] <- t(cor_mat)[lower.tri(cor_mat)]
  
  # Set the diagonal of the correlation matrix to 1
  diag(cor_mat) <- 1
  
  # Compute the standard error matrix using the diagonal variance values
  se_mat <- diag(sqrt(exp(v_vec[diag_index])))
  
  # Compute the covariance matrix by combining the correlation and standard error matrices
  V_mat <- se_mat %*% cor_mat %*% se_mat
  
  return(V_mat)
}


#' @title The function computes the EB1 and EB2 matrices based on the provided input parameters.
#'
#' @param alpha Input alpha matrix.
#' @param mu1 Input mu1 list.
#' @param mu2 Input mu2 list.
#' @param column_config Configuration of columns.
#' @param nsex Number of sexes.
#' @param nsnp Number of SNPs.
#'
#' @return A list containing the computed EB1 and EB2 matrices.
#' @export

compute_b1b2<-function(alpha,mu1,mu2,column_config,nsex,nsnp){

  # Initialize matrices with zeros
  EB1 = matrix(0, ncol = nsex, nrow = nsnp)
  EB2 = matrix(0, ncol = 3, nrow = nsnp)
   
  # Loop over column configurations
  for(x in 1:length(column_config)){
	cor_cols = column_config[[x]]
		
	# Update EB1 and EB2 matrices based on column configurations
	EB1[,cor_cols] = EB1[,cor_cols] + alpha[,x] * mu1[[x]]		
	if(length(cor_cols) == 2){
	  EB2 = EB2 + alpha[,x] * mu2[[x]]
	}else{
	  EB2[,cor_cols] = EB2[,cor_cols] + alpha[,x] * mu2[[x]]
	}
  }

	return(list(EB1 = EB1, EB2 = EB2))

}


#' @title The function computes the posterior inclusion probability for the l'th effect using log Bayes Factors.
#'
#' @param lbf Log Bayes Factors.
#' @param prior_weights A vector of prior weights.
#'
#' @return A list containing the posterior inclusion probability (alpha_wmulti) and the logarithm of the likelihood (loglik).
#' @export

compute_softmax <- function(lbf, prior_weights) {
  
  # Subtracting max lbf from all lbf values for numerical stability
  maxlbf = max(lbf)
  w_multi = exp(lbf - maxlbf)
  
  # Multiply by prior weights
  w_weighted_multi = w_multi * prior_weights
  
  # Sum over weighted values
  weighted_sum_w = sum(w_weighted_multi)
  
  # Compute the normalized weights
  alpha_wmulti = w_weighted_multi / weighted_sum_w
  
  # Return computed values
  return(list(alpha_wmulti = alpha_wmulti, loglik = maxlbf + log(weighted_sum_w)))
}


#' @title The function calculates the number of SNPs in a credible set required to achieve the desired coverage.
#'
#' @param x A numeric vector representing SNP values.
#' @param coverage The desired coverage (default is 0.95 or 95%).
#'
#' @return The number of SNPs needed in the credible set.
#' @export

n_in_CS_x <- function (x, coverage = 0.95) {
  
  # Cumulatively sum the sorted SNPs in descending order 
  # and find the position just before the cumulative sum exceeds the desired coverage
  num_snps_required = sum(cumsum(sort(x, decreasing = TRUE)) < coverage) + 1
  
  return(num_snps_required)
}


#' @title The function identifies the SNPs that are included in the credible set to achieve the desired coverage.
#'
#' @param x A numeric vector representing SNP values.
#' @param coverage The desired coverage (default is 0.95 or 95%).
#'
#' @return A binary vector, where 1 indicates the SNP is in the credible set and 0 otherwise.
#' @export

in_CS_x <- function (x, coverage = 0.95) {
  
  # Get the number of SNPs required for the desired coverage
  n = n_in_CS_x(x, coverage)
  
  # Order the SNPs in descending order
  o = order(x, decreasing = TRUE)
  
  # Initialize result vector as zeros
  result = rep(0, length(x))
  
  # Mark the top 'n' SNPs as 1, indicating they are in the credible set
  result[o[1:n]] = 1
  
  return(result)
}


#' @title The function applies the in_CS_x function row-wise to a matrix to identify the SNPs that are included in the credible set for each row.
#'
#' @param res A matrix where each row represents a set of SNP values.
#' @param coverage The desired coverage (default is 0.95 or 95%).
#'
#' @return A matrix where each row is a binary vector indicating which SNPs are in the credible set.
#' @export

in_CS <- function (res, coverage = 0.95) {
  
  # Apply the in_CS_x function to each row of the 'res' matrix
  result_matrix = t(apply(res, 1, function(x) in_CS_x(x, coverage)))
  
  return(result_matrix)
}


#' @title The function computes the minimum, mean, and median of absolute correlations for specific positions.
#'
#' @param pos A vector of positions to consider.
#' @param Xcorr A correlation matrix.
#'
#' @return A vector containing the minimum, mean, and median of the absolute correlations.
#' @export

get_purity <- function(pos, Xcorr) {
  
  # If there's only one position, return a vector of ones
  if (length(pos) == 1) {
    return(c(1, 1, 1))
  } else {
    # Extract absolute correlations for specified positions
    value_list = lapply(Xcorr, function(x) c(abs(x[pos, pos])))
    
    # Combine the lists column-wise
    value_matrix = Reduce(cbind, value_list)
    
    # Compute the maximum for each row
    value_max = do.call(pmax, data.frame(value_matrix))
    
    # Return the desired statistics
    return(c(min(value_max, na.rm = TRUE),
             mean(value_max, na.rm = TRUE),
             median(value_max, na.rm = TRUE)))
  }
}


#' @title The function extracts credible sets from sdSuSiE results.
#'
#' @param res The result object from sdSuSiE.
#' @param Xcorr A list of correlation matrix used in sdSuSiE.
#' @param coverage A numeric value indicating the desired coverage for the credible set, default is \(0.95\).
#' @param prior_tol A numeric threshold value to determine prior variances, default is \(1e-9\).
#' @param cor_method A string indicating which method to use when computing purity, defalt is "min.abs.corr".
#' @param cor_threshold A numeric threshold value for purity measurement, default is 0.5.
#'
#' @return A list containing:
#' \item{cs}{The credible sets.}
#' \item{cs_category}{The category of each credible set having either common effect, sex-dimorphic effect, or both.}
#' \item{purity}{Purity statistics for each credible set.}
#' \item{cs_index}{The index for the credible sets.}
#' \item{coverage}{The achieved coverage for the credible sets.}
#' \item{requested_coverage}{The desired coverage (provided as input).}
#' 
#' @export

sdSuSiE_get_cs <- function(res, Xcorr, coverage = 0.95, prior_tol = 1e-9, cor_method = cor_method, cor_threshold = cor_threshold){
  
  include_idx = unlist(lapply(res$V, function(x)max(diag(x)) > prior_tol))
  alpha_either = t(Reduce(cbind, lapply(res$alpha, function(x)apply(x, 1, sum))))
  status = in_CS(alpha_either, coverage = 0.95)
  cs = lapply(1:nrow(status), function(i) which(status[i,] != 0))
  claimed_coverage = sapply(1:length(cs), function (i) sum(alpha_either[i,][cs[[i]]]))
  include_idx = include_idx * (lapply(cs, length) > 0)
  include_idx = include_idx * (!(duplicated(cs)))
  
  include_idx = as.logical(include_idx)
  if (sum(include_idx) == 0)
    return(list(cs = NULL, cs_category = NULL, purity = NULL, cs_index = NULL, coverage = NULL, requested_coverage = coverage))
  
  cs = cs[include_idx]
  claimed_coverage = claimed_coverage[include_idx]
  
  purity = data.frame(do.call(rbind, lapply(1:length(cs), function (i){
    get_purity(cs[[i]], Xcorr)
  })))
  
  colnames(purity) = c("min.abs.corr", "mean.abs.corr", "median.abs.corr")
  threshold = cor_threshold
  is_pure = which(purity[,colnames(purity) %in% cor_method] >= threshold)
  
  if (length(is_pure) > 0) {
    cs = cs[is_pure]
    purity = purity[is_pure,]
    claimed_coverage = claimed_coverage[is_pure]
    row_names = paste0("L",which(include_idx)[is_pure])
    names(cs) = row_names
    rownames(purity) = row_names
    
    cs_index = which(include_idx)[is_pure]
    cs_category = unlist(lapply(res$alpha[cs_index], function(x)res$name_config[which.max(colSums(x))]))
    names(cs_category) = row_names
    
    return(list(cs = cs,
                cs_category = cs_category,
                purity   = purity,
                cs_index = cs_index,
                coverage = claimed_coverage,
                requested_coverage=coverage))
  } else
    return(list(cs = NULL, cs_category = NULL, purity = NULL, cs_index = NULL, coverage = NULL, requested_coverage = coverage))
}


#' @title The function extracts the maximum posterior inclusion probabilities (PIP) across all components from the given results.
#'
#' @param res The result object from sdSuSiE.
#'
#' @return A vector containing the maximum PIP for each SNP.
#' 
#' @export

sdSuSiE_get_pip_either <- function(res) {
  
  # Aggregate the alpha values for each subset
  aggregated_alpha = Reduce(cbind, lapply(res$alpha, function(x) apply(x, 1, sum)))
  
  # Return the maximum value for each result
  return(do.call(pmax, data.frame(aggregated_alpha)))
}


#' @title The function extracts the maximum posterior inclusion probabilities (PIP) for each configuration from the given results.
#'
#' @param res The result object from sdSuSiE.
#'
#' @return A matrix containing the maximum PIP for each causal configuration.
#' 
#' @export

sdSuSiE_get_pip_config <- function(res) {
  
  # Extract the maximum alpha value for each column across all alpha matrices
  max_alpha_per_config = Reduce(cbind, lapply(1:ncol(res$alpha[[1]]), function(x) {
    do.call(pmax, lapply(res$alpha, function(y) y[,x]))
  }))
  
  return(max_alpha_per_config)
}


#' @title The function creates GWAS plots for sex-stratified GWASs.
#'
#' @param data_plot A data frame containing results from sex-stratified GWASs.
#'
#' @return A ggplot object visualizing the results of sex-stratified GWASs. 

#' @import dplyr ggplot2 ggrepel 
#' @export

gwas_plot_fundouble <- function(data_plot) { 
  p_manhattan <- ggplot() +
	geom_point(data = data_plot %>% filter(Index_SNP == 0 & sex == "M"), aes(x = POS, y = P, color = r2), size = 2) +
    geom_point(data = data_plot %>% filter(Index_SNP == 1 & sex == "M"), aes(x = POS, y = P), size = 3, shape = 18, color = "red") +
    geom_text_repel(data = data_plot %>% filter(Index_SNP == 1 & sex == "M"), aes(x = POS, y = P, label = SNP), vjust = 1.2, size = 4, show.legend = FALSE) +
	geom_point(data = data_plot %>% filter(Index_SNP == 0 & sex == "F"), aes(x = POS, y = P, color = r2), size = 2) +
    geom_point(data = data_plot %>% filter(Index_SNP == 1 & sex == "F"), aes(x = POS, y = P), size = 3, shape = 18, color = "red") +
    geom_text_repel(data = data_plot %>% filter(Index_SNP == 1 & sex == "F"), aes(x = POS, y = P, label = SNP), vjust = 1.2, size = 4, show.legend = FALSE) +
    scale_color_stepsn(
      colors = c("navy", "lightskyblue", "green", "orange", "red"),
      breaks = seq(0.2, 0.8, by = 0.2),
      limits = c(0, 1),
      show.limits = TRUE,
      na.value = 'grey50',
      name = expression(R^2)
    ) +
    geom_hline(
      yintercept = c(0,-log10(5e-8),log10(5e-8)),
      linetype = "dashed",
      color = "grey50",
      size = 0.5
    ) +
    geom_vline(
      xintercept = data_plot %>% filter(Index_SNP == 1) %>% pull(POS),
      linetype = "dashed",
      color = "grey50",
      size = 0.5
    ) +
    xlim(min(data_plot$POS), max(data_plot$POS)) +
    #expand_limits(x = round(max(data_plot$POS)/1e6)*1e6) +
    scale_x_continuous(labels = function(x) {
      if (max(x, na.rm = TRUE) > 1e6) {
        paste0(x / 1e6, " MB")
      } else {
        paste0(x / 1e3, " KB")
      }
    }) +
	labs(x = "", y = "-log10(p-value)", title = "Male vs Female GWAS") +
    guides(fill = guide_legend(title = as.expression(bquote(R^2)))) +
    theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(color = "black", size = 10, hjust = 0.5),
    panel.grid = element_blank(),
    axis.line.y = element_line(color = "black", linetype = "solid"),
    axis.line.x = element_line(color = "black", linetype = "solid"),
    panel.border = element_rect(linetype = "solid", fill = NA),
    panel.background = element_blank(),
    legend.position = c(0.1, 0.78),
    legend.box.background = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(color = "black", size = 8),
	legend.title=element_text(color = "black", size = 8)
  )
  
  return(p_manhattan)
}


#' @title The function creates a GWAS or PIP plot for visualizing results from either univariate sex-dimorphic analysis or sdSuSiE.
#'
#' @param data_plot A data frame containing results from sex-stratified GWASs or sdSuSiE.
#' @param title_name A string specifying the title of the plot.
#' @param ylab_name A string specifying the label for the y-axis, such as "-log10(p-value)" or "PIP".
#' @param yintercept A numeric value indicating the significance line, such as the genome-wide significance threshold.
#'
#' @return A ggplot object visualizing the results of either univariate sex-dimorphic analysis or sdSuSiE.
#'
#' @import dplyr ggplot2 ggrepel
#' @export

gwas_plot_fun <- function(data_plot, title_name, ylab_name, yintercept) { 
  
  p_manhattan <- ggplot() +
	geom_point(data = data_plot %>% filter(Index_SNP == 0), aes(x = POS, y = P_PIP, color = r2), size = 2) +
    geom_point(data = data_plot %>% filter(Index_SNP == 1), aes(x = POS, y = P_PIP), size = 3, shape = 18, color = "red") +
    geom_text_repel(data = data_plot %>% filter(Index_SNP == 1), aes(x = POS, y = P_PIP, label = SNP), vjust = 1.2, size = 4, show.legend = FALSE) +
    scale_color_stepsn(
      colors = c("navy", "lightskyblue", "green", "orange", "red"),
      breaks = seq(0.2, 0.8, by = 0.2),
      limits = c(0, 1),
      show.limits = TRUE,
      na.value = 'grey50',
      name = expression(R^2)
    ) +
    geom_hline(
      yintercept = yintercept,
      linetype = "dashed",
      color = "grey50",
      size = 0.5
    ) +
    geom_vline(
      xintercept = data_plot %>% filter(Index_SNP == 1) %>% pull(POS),
      linetype = "dashed",
      color = "grey50",
      size = 0.5
    ) +
    xlim(min(data_plot$POS), max(data_plot$POS)) +
    #expand_limits(x = round(max(data_plot$POS)/1e6)*1e6) +
    scale_x_continuous(labels = function(x) {
      if (max(x, na.rm = TRUE) > 1e6) {
        paste0(x / 1e6, " MB")
      } else {
        paste0(x / 1e3, " KB")
      }
    }) +
	labs(x = "", y = ylab_name, title = title_name) +
    guides(fill = guide_legend(title = as.expression(bquote(R^2)))) +
    theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(color = "black", size = 10, hjust = 0.5),
    panel.grid = element_blank(),
    axis.line.y = element_line(color = "black", linetype = "solid"),
    axis.line.x = element_line(color = "black", linetype = "solid"),
    panel.border = element_rect(linetype = "solid", fill = NA),
    panel.background = element_blank(),
    legend.position = "none",
    legend.box.background = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(color = "black", size = 8)
  )
  
  return(p_manhattan)
}


#' @title The function creates GWAS plots from sex-stratified GWASs and univariate sex-dimorphic analysis, and PIP plots from sdSuSiE.
#'
#' @param res The result object from sdSuSiE.
#' @param R_mat A list of length 2 containing LD matrices for each sex.
#' @param summary_data A list of data frames for each sex containing at least SNP, Z, and possibly POS (position) columns.
#'
#' @return A combined ggplot object visualizing the results of sex-stratified GWASs, univariate sex-dimorphic analysis and sdSuSiE.
#'
#' @import ggplot2 ggrepel cowplot
#' @export

sdSuSiE_Plot <- function(res, R_mat, summary_data) {
  
  # Extract names of the datasets in summary_data
  name_vec <- names(summary_data)
  
  # Check if Position information is provided; if not, generate a position based on row order
  if (!("POS" %in% colnames(summary_data[[name_vec[1]]]))) {
    summary_data <- lapply(summary_data, function(x) {
      x$POS = seq_len(nrow(x))
      return(x)
    })
  }
  
  # plot_a
  # Identify lead SNP in each sex
  target_rs = c(summary_data[[1]]$SNP[which.max(abs(summary_data[[1]]$Z))],summary_data[[2]]$SNP[which.max(abs(summary_data[[2]]$Z))])
  # Extract LD (r-squared) 
  r2_M = R_mat[[1]][,colnames(R_mat[[1]])==target_rs[1]]^2
  r2_F = R_mat[[2]][,colnames(R_mat[[2]])==target_rs[2]]^2
  # Construct the data frame for GWAS plots based on the sex-stratified GWASs
  SNP = summary_data[[1]]$SNP
  POS = as.numeric(summary_data[[1]]$POS)
  P_M = -log10(2*pnorm(-abs(summary_data[[1]]$Z)))
  P_F = -log10(2*pnorm(-abs(summary_data[[2]]$Z)))
  index_SNP_col = rep(0,length(SNP))
  index_SNP_col[which(SNP %in% target_rs[1])] = 1
  GWAS_plot_M <- data.frame(SNP, POS, r2_M, P_M, index_SNP_col)
  colnames(GWAS_plot_M)[3:5]=c("r2","P","Index_SNP")
  GWAS_plot_M$sex="M" 
  index_SNP_col = rep(0,length(SNP))
  index_SNP_col[which(SNP %in% target_rs[2])] = 1
  GWAS_plot_F <- data.frame(SNP, POS, r2_F, P_F, index_SNP_col)
  colnames(GWAS_plot_F)[3:5] = c("r2","P","Index_SNP")
  GWAS_plot_F$sex = "F"
  GWAS_plot_F$P = -1*P_F
  GWAS_plot = rbind(GWAS_plot_M, GWAS_plot_F)
  plot_a = gwas_plot_fundouble(GWAS_plot)

  # plot_b
  # Calculate Z-score in the univariate sex-dimorphic analysis
  GWAS_sd = GWAS_plot_M
  GWAS_sd$Z = (summary_data[[1]]$Beta - summary_data[[2]]$Beta)/sqrt(summary_data[[1]]$Se^2 + summary_data[[2]]$Se^2)
  # Identify lead SNP in each sex
  lead_SNP = GWAS_sd$SNP[which.max(abs(GWAS_sd$Z))]
  # Extract LD (r-squared)
  r2 = ((summary_data[[1]]$N*R_mat[[1]] + summary_data[[2]]$N*R_mat[[2]])/(summary_data[[1]]$N + summary_data[[2]]$N))[,which(summary_data[[1]]$SNP==lead_SNP)]^2
  P_PIP = -log10(2*pnorm(-abs(GWAS_sd$Z)))
  lead_SNP_col=rep(0,length(SNP))
  lead_SNP_col[which(SNP %in% lead_SNP)]=1
  GWAS_sd_plot <- data.frame(SNP, POS, r2, P_PIP, lead_SNP_col)
  colnames(GWAS_sd_plot)[5]="Index_SNP"
  plot_b <- gwas_plot_fun(GWAS_sd_plot, "Univariate Sex-dimorphic Analysis", "-log10(p-value)", -log10(5e-8))

  # plot_c
  # PIP_c from sdSuSiE
  GWAS_sd_plot$P_PIP=res$PIP_config[,1]
  GWAS_sd_plot$Index_SNP = 0
  GWAS_sd_plot$Index_SNP[which.max(res$PIP_config[,1])] = 1
  plot_c <- gwas_plot_fun(GWAS_sd_plot, "PIPc from sdSuSiE", "PIP", 0.5)
  
  # plot_d
  # PIP_d from sdSuSiE
  GWAS_sd_plot$P_PIP=res$PIP_config[,2]
  GWAS_sd_plot$Index_SNP = 0
  GWAS_sd_plot$Index_SNP[which.max(res$PIP_config[,2])] = 1
  plot_d <- gwas_plot_fun(GWAS_sd_plot, "PIPd from sdSuSiE", "PIP", 0.5)
  
  # plot_e
  # PIP_b from sdSuSiE
  GWAS_sd_plot$P_PIP=res$PIP_config[,3]
  GWAS_sd_plot$Index_SNP = 0
  GWAS_sd_plot$Index_SNP[which.max(res$PIP_config[,3])] = 1
  plot_e <- gwas_plot_fun(GWAS_sd_plot, "PIPb from sdSuSiE", "PIP", 0.5)

  # Combine all plots and return
  col1 <- plot_grid(plot_a, plot_b, labels = c("a", "b"), hjust = -1, vjust = 1.5, nrow = 2,ncol = 1, rel_heights = c(2,1))
  col2 <- plot_grid(plot_c, plot_d, plot_e, labels = c("c", "d", "e"), hjust = -1, vjust = 1.5, nrow = 3, ncol = 1)
  combined_plot <- plot_grid(col1, col2, ncol = 2, nrow = 1)
  
  return(combined_plot)
}
