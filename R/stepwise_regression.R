#' @title The function computes the inverse of a matrix.
#'
#' @description This function attempts to compute the inverse of a matrix by first using Cholesky decomposition. 
#' If Cholesky decomposition fails, it falls back to LU decomposition. If both methods fail, it returns an error message.
#'
#' @param V A square, symmetric, positive-definite matrix for which the inverse is to be computed.
#'
#' @return A list with the following elements:
#' \item{success}{A logical value indicating whether the inverse computation was successful.}
#' \item{inverse}{The computed inverse matrix, if successful. Returned as a matrix.}
#' \item{error}{An error message if both decompositions fail. Returned as a character string.}
#'
#' @export

comput_inverseR <- function(V) {

  tryCatch({
    Vc <- chol(V) 
    invVc <- solve(Vc) 
    invV <- invVc %*% t(invVc) 
    return(list(success = TRUE, inverse = invV))
  }, error = function(e) {
    message("Cholesky decomposition failed. Attempting LU decomposition.")
    
  tryCatch({
	lu_decomp <- qr(V)
	L <- qr.R(lu_decomp) 
	P <- qr.Q(lu_decomp) 
      
	invU <- solve(L) 
	invL <- solve(P) 
	invV <- invU %*% invL 
	return(list(success = TRUE, inverse = invV))
  }, error = function(e) {
      message("LU decomposition failed. The matrix might be singular.")
      return(list(success = FALSE, error = "Both Cholesky and LU decomposition failed"))
    })
  })
}


#' @title The main function for stepwise regression in sex-dimorphic mapping.
#'
#' @description The stepwise regression extends the standard stepwise selection procedure in the COJO framework to fine-map sex-dimorphic SNPs. 
#'
#' @param summa A data frame containing summary statistics of SNPs for marginal common and sex-dimorphic effects, including columns such as `SNP`, `Beta`, `Se`, `Z`, `N`, and `p`.
#' These summary statistics can be derived using the sex-stratified GWAS summary statistics of males and females. The names of SNP with marginal effects are prefixed with "S*".
#' @param LD A 2p * 2p matrix. It represents pairwise correlations among SNPs with both common and sex-dimorphic effects.
#' @param p_cutoff A numeric value specifying the p-value threshold for selecting significant SNPs in the stepwise regression process. Default is 5e-8.
#' @param collinear A numeric value specifying the maximum allowable collinearity (1 - 1/diag(inv_sel)) among selected SNPs. Default is 0.9.
#' @param ncore An integer specifying the number of cores to use for parallel processing. Default is 10.
#' 
#' @return A data frame containing information about the selected SNPs, including the following columns:
#' \item{b_c}{The effect size estimates of selected SNPs in the jointly fitted model.}
#' \item{se_c}{The standard errors of the effect size estimates.}
#' \item{p_c}{The p-values of selected SNPs in the jointly fitted model.}
#' \item{LD_r}{The pairwise LD values among selected SNPs.}
#' \item{sig}{Indicator of significance for each SNP (1 if `p_c < 5e-8`, 0 otherwise).}
#'
#' If no significant SNPs are identified, the function returns `NULL`.
#'
#' @import parallel doParallel foreach
#' @export
#' @examples
#' library(sdSuSiE)
#' data(GWAS_M)
#' data(GWAS_F)
#' data(LD_M)
#' data(LD_F)
#' SNP = GWAS_M$SNP
#' Z = (sqrt(GWAS_M$N - 1) * GWAS_M$Z + sqrt(GWAS_F$N - 1) * GWAS_F$Z) / sqrt(GWAS_M$N + GWAS_F$N - 1)
#' N = round(GWAS_M$N + GWAS_F$N)
#' Beta = Z / sqrt(N - 1)
#' Se = 1 / sqrt(N - 1)
#' data1 = data.frame(SNP, Beta, Se, Z, N)
#' SNP = paste0("S*", SNP)
#' Z = (-sqrt(GWAS_M$N - 1) * GWAS_M$Z + sqrt(GWAS_F$N - 1) * GWAS_F$Z) / sqrt(GWAS_M$N + GWAS_F$N - 1)
#' N = round(GWAS_M$N + GWAS_F$N)
#' Beta = Z / sqrt(N - 1)
#' Se = 1 / sqrt(N - 1)
#' data2 = data.frame(SNP, Beta, Se, Z, N)
#' data = rbind(data1, data2)
#' data$p = pchisq((data$Beta/data$Se)^2, df = 1, lower.tail = FALSE)
#' N_M = round(median(GWAS_M$N))
#' N_F = round(median(GWAS_F$N))
#' N = N_M + N_F
#' V1 = (N_M-1)*LD_M + (N_F-1)*LD_F
#' V2 = (-1)*(N_M-1)*LD_M + (N_F-1)*LD_F
#' LD = rbind(cbind(V1,V2),cbind(V2,V1)) / (N-1)
#' res = sdstepwise(data, LD, ncore = 1)

sdstepwise <- function(summa, LD, p_cutoff = 5e-8, collinear = 0.9, ncore = 10){
  is_run=TRUE
  ##adjustment for effective sample size
  #h=2*summa$freq*(1-summa$freq)
  #vary=median(h*summa$Beta^2+h*summa$Se^2*summa$N)
  #summa$n=(vary - h*summa$Beta^2)/(h*summa$Se^2) + 1.0
  #summa$Z=summa$Beta/summa$Se
  #summa$Se=sqrt(1/(summa$Z^2+summa$n))
  #summa$Beta=summa$Z*summa$Se
  
  #Obtain the jointly model p-value for each SNP
  cl <- makeCluster(min(ncore, detectCores() - 1))
  registerDoParallel(cl)
  results <- foreach(i = 1:(dim(LD)[1]/2), .combine = 'rbind', .export = c('comput_inverseR')) %dopar% {
	tidx = c(match(summa$SNP[i], summa$SNP), match(summa$SNP[i], summa$SNP) + (dim(LD)[1]/2))            
	inv_sel = comput_inverseR(as.matrix(LD[tidx, tidx]))$inverse 
	b_cond = inv_sel %*% summa$Beta[tidx]                     
	se_cond = sqrt(diag(inv_sel) / summa$N[tidx])          
	p_cond = pchisq((b_cond / se_cond)^2, df = 1, lower.tail = FALSE) 
	c(p_c = p_cond[1], p_d = p_cond[2])                    
  }
  stopCluster(cl)
  summa$p_i = c(results[, 1],results[, 2])
  summa$idx = c(1:length(summa$p))
  summa$SNPu = gsub("S\\*","",summa$SNP)
  
  #Select the SNP with the lowest p-value in the jointly model for each SNP
  select = grep(summa$SNPu[which.min(summa$p_i)], summa$SNPu)

  #Check whether SNPs with significant effects, if not return to NULL
  if(min(summa$p_i) < p_cutoff){
    #Stepwise selection
    while(is_run){
	  #Store the last selected variables 
	  selectt = select
	  unselect = c(1:dim(summa)[1])[!c(1:dim(summa)[1]) %in% select]
	  inv_sel = comput_inverseR(as.matrix(LD[select,select]))$inverse
	  is_run = FALSE

	  #Check the collinear issue
	  if(max(1 - 1/diag(inv_sel)) > collinear){
	    tmp = summa[!summa$idx %in% select,]
	    if(min(tmp$p_i) < p_cutoff){
		  select[c(length(select)-1, length(select))] = grep(tmp$SNPu[which.min(tmp$p_i)], tmp$SNPu)
	    }
	  }else{
	    #If there no collinear issue, perform the conditional analysis
	    cl <- makeCluster(ncore) 
	    registerDoParallel(cl)
	    p_value <- foreach(i = unselect, .combine = c, .export = c('comput_inverseR')) %dopar% {
	      tidx = c(select, grep(summa$SNPu[i], summa$SNPu))
	      tidx = tidx[tidx != i]
 
		  inv_sel = comput_inverseR(as.matrix(LD[tidx, tidx]))$inverse
		
		  if(!is.null(inv_sel)){
		    RCR = inv_sel %*% LD[tidx, i]
  
		    b_value = summa$Beta[i] - sum(RCR * summa$Beta[tidx])
		    se_value = sqrt((1 - sum(LD[tidx, i] * RCR)) / summa$N[i])
		    pchisq((b_value / se_value)^2, df = 1, lower.tail = FALSE)
		  }else{
		    NA
		  }
	    }
	    stopCluster(cl)

	    #If there exsits the conditional p-value below than the p_cutoff, add the SNP to the select SNPs
	    is_p = TRUE
  	    tmp = data.frame(unselect, p_value)
	    while(is_p){
	      is_p = FALSE
	      if(min(na.omit(tmp$p_value)) < p_cutoff){
		    #Check the collinear issue, if there exists the collinear issue, drop the new added SNP, and select the next largest one if its below than the p_cutoff
		    ps = grep(summa$SNPu[tmp$unselect[which.min(tmp$p_value)]], summa$SNPu)
		    m = comput_inverseR(as.matrix(LD[c(select, ps), c(select, ps)]))$inverse
		    if(max(1-1/diag(m)) <= collinear & !is.null(m)){
		      select = c(select, ps)
		    }else{
		      de = which.min(tmp$p_value)
		      if(de <= dim(summa)[1]/2){
			    de = c(de, de + dim(summa)[1]/2)
		      }else{
			    de=c(de - dim(summa)[1]/2, de)
		      }
		      tmp = tmp[-de,]
		      is_p = TRUE
		    }
	      }
	    }
	  }

	  #Fit all select SNPs jointly in the model and drop the SNP with largest p-value that is greater than the p_cutoff
	  if(length(select) > length(selectt)){
	    inv_sel = comput_inverseR(as.matrix(LD[select, select]))$inverse
	    b_c = inv_sel %*% summa$Beta[select]
	    se_c = sqrt(diag(inv_sel)/summa$N[select])
	    p_c <- pchisq((b_c/se_c)^2, df = 1, lower.tail = FALSE)
	    p_n <- sapply(seq(1, length(p_c), by = 2), function(i) min(p_c[i], p_c[i + 1]))
	    if(max(p_n) > p_cutoff){
	      select = select[-c(2*which.max(p_n)-1, 2*which.max(p_n))]
	    }
	  }
	  if(!setequal(select, selectt)){
	    is_run=TRUE
	  }
    }
  
    #Perform the jointly analysis
    inv_sel = comput_inverseR(as.matrix(LD[select, select]))$inverse
    b_c = inv_sel %*% summa$Beta[select]
    se_c = sqrt(diag(inv_sel)/summa$N[select])
    p_c <- pchisq((b_c/se_c)^2, df = 1, lower.tail = FALSE)
    LD_r = rep(0, length(select))
    if(length(select) > 1){
	  for(i in 1:(length(select)-1)){
	    LD_r[i] = LD[select[i], select[i+1]]
	  }
    }
    data = data.frame(summa[select,], b_c, se_c, p_c, LD_r)
    data$sig = 0
    data$sig[data$p_c < 5e-8] = 1
    return(data)
  }
  else{
	return(NULL)
  }
}
