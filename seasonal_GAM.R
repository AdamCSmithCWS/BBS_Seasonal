
### a function to calculate the basis function(s) for a GAM component of a JAGS model
### orig.preds = predictor value on which the GAM smooth is modeled, e.g., a column in the data.frame with 1-row for each statistical unit
### nknots = number of knots in the GAM (careful! you need to know what you're doing :-) 
### random = logical; TRUE if the GAM parameters are to be estimated as random effects (e.g., separate smooths within strata, centered on a mean smooth across all strata)
### predpoints = optional vector of points along the range of hte predictor at which the user wants to generate predictions
### standardize = logical; TRUE if the predictor should be standardized (centered and re-scaled)






gam.basis.func <- function(orig.preds = dd[,"julian.date"],
                           nknots = 12,
                           standardize = T,
                           random = T,
                           predpoints = NULL){
  
  if(any(is.na(orig.preds) == T)){
    stop("This GAM formulation cannot handle missing values in the predictor")
  }
  
  
  vmin = min(orig.preds) 
  vmax = max(orig.preds)
  
  
  
  
  if(is.null(predpoints) == F){
    npredpoints = length(predpoints) 
  }else{
    npredpoints = NA
  }
  
  if(standardize){ 
    
    recenter = floor(diff(c(vmin,vmax))/2) #centering to the middle of the range
    #ignores the distribution of the data points, which may be dangerous of the data are not well distributed
    rescale = ((vmax-vmin)+1)/2 # rescales to a unit value = 1/2 of the range, similar to SD = 1 for BBS data
    vscale = (orig.preds-recenter)/rescale
  }else{
    vscale = orig.preds
  }
  
  
  
  yminsc = min(vscale,na.rm = T)
  ymaxsc = max(vscale,na.rm = T)
  
  ############ GAM basis function
  
  knotsgamx<- seq(yminsc,ymaxsc,length=(nknots+2))[-c(1,nknots+2)]
  
  gamx_OMEGA_all<-(abs(outer(knotsgamx,knotsgamx,"-")))^3
  gamx_svd.OMEGA_all<-svd(gamx_OMEGA_all)
  gamx_sqrt.OMEGA_all<-t(gamx_svd.OMEGA_all$v  %*% (t(gamx_svd.OMEGA_all$u)*sqrt(gamx_svd.OMEGA_all$d)))
  
  if(is.null(predpoints) == F){
    
    gamx_K<-(abs(outer(seq(yminsc,ymaxsc,length = npredpoints),knotsgamx,"-")))^3
    gamx.basispred<-t(solve(gamx_sqrt.OMEGA_all,t(gamx_K)))
  }else{
    gamx.basispred = NA
  }
  
  gamx_K1<-(abs(outer(vscale,knotsgamx,"-")))^3
  gamx.basis<-t(solve(gamx_sqrt.OMEGA_all,t(gamx_K1)))
  
  if(random){
    mod.code <- "
    ###########insert the following text into an existing JAGS model##############
    ###########making sure to link the relevant parameters into section##############
    ###########of the model that estimates the full likelihood ##############
    ###################################################################
    # the GAM smooth is predicted as intercepts for each whole value of the
    # predictor (e.g., each year in a temporal smooth)
    # so the likelihood line will need to include a component such as the following for each data point-i
    # y[i] <- ..... + gam.sm[i] +....
    
    for ( i in 1:n ){
    for(k in 1:nknots){
    gamx.part[i,k] <- beta.gamx[group[i],k]*(gamx.basis[i,k])
    }
    gam.sm[i] <- sum(gamx.part[i,1:nknots])
    }#i
    
    
    ###################################################################
    ##### this text assumes the GAM parameters (betas) are estimated as random effects########
    ##### your will need to provide the variable ngroups (an integer indicating the number of groups in the random effect)####
    ##### and ensure that your R-script retains the correct order of groups in relation to the predicted values indexed by the 1:ngroups section below ######
    
    # Following Crainiceanu, C. M., Ruppert, D. & Wand, M. P. (2005). Bayesian Analysis for Penalized Spline Regression Using WinBUGS. Journal of Statistical Softare, 14 (14), 1-24.
    
    taugamx~dgamma(1.0E-2,1.0E-4) #alternate prior, original from Cainiceanu et al. second gamma parameter == 0.0001 << (abs(mean(B.gamx[]))^2)/2, mean(B.gamx[]) ~ 0.2
    sdgamx <- 1/(pow(taugamx,0.5)) # 
    taubeta <- 1/pow(sdbeta,2) # prior on precision of gam coefficients(
    sdbeta ~ dunif(0,5)
    
    for(k in 1:nknots)
    { # Computation of GAM components
    B.gamx[k] ~ dnorm(0,taugamx)
    
    for(j in 1:ngroups)
    {
    beta.gamx[j,k] ~ dnorm(B.gamx[k],taubeta)
    
    
    for ( i in 1:npredpoints ){ #this allows the user to model just the GAM component at particular values of the predictor
    gamx.partpred[i,k,j] <- beta.gamx[j,k]*(gamx.basispred[i,k])
    
    }#i
    
    }#k
    }#j
    
    for (i in 1:npredpoints){
    for(j in 1:ngroups)
    {
    gam.smooth[i,j] <- sum(gamx.part[i,1:nknots,j])
    }#k
    }#i
    "
  }else{
    mod.code = "
    ###########insert the following text into an existing JAGS model##############
    ###########making sure to link the relevant parameters into section##############
    ###########of the model that estimates the full likelihood ##############
    ###################################################################
    ###################################################################
    # the GAM smooth is predicted as intercepts for each whole value of the
    # predictor (e.g., each year in a temporal smooth)
    # so the likelihood line will need to include a component such as the following for each data point-i
    # y[i] <- ..... + gam.sm[i] +....
    
    for ( i in 1:n ){
    for(k in 1:nknots){
    gamx.part[i,k] <- beta.gamx[k]*(gamx.basis[i,k])
    }
    gam.sm[i] <- sum(gamx.part[i,1:nknots])
    }#i
    
    
    ###################################################################
    ##### this text assumes the GAM parameters (betas) are estimated as random effects########
    ##### your will need to provide the variable ngroups (an integer indicating the number of groups in the random effect)####
    ##### and ensure that your R-script retains the correct order of groups in relation to the predicted values indexed by the 1:ngroups section below ######
    
    # Following Crainiceanu, C. M., Ruppert, D. & Wand, M. P. (2005). Bayesian Analysis for Penalized Spline Regression Using WinBUGS. Journal of Statistical Softare, 14 (14), 1-24.
    
    taugamx~dgamma(1.0E-2,1.0E-4) #alternate prior, original from Cainiceanu et al. second gamma parameter == 0.0001 << (abs(mean(B.gamx[]))^2)/2, mean(B.gamx[]) ~ 0.2
    sdgamx <- 1/(pow(taugamx,0.5)) # 
    
    for(k in 1:nknots)
    { # Computation of GAM components
    beta.gamx[k] ~ dnorm(0,taugamx)
    
    
    for ( i in 1:npredpoints )
    {
    gamx.partpred[i,k] <- beta.gamx[k]*(gamx.basispred[i,k])
    
    }#i
    
    }#k
    
    
    for (i in 1:npredpoints)
    {
    
    gam.smooth[i] <- sum(gamx.part[i,1:nknots])
  }#i
    " 
}
  
  writeLines(mod.code,con = "text to add to jags model.txt")
  
  return(list(gamx.basis = gamx.basis,
              gamx.basispred = gamx.basispred,
              vscale = vscale,
              orig.preds = orig.preds,
              predpoints = predpoints,
              nknots = nknots,
              npredpoints = npredpoints,
              mod.code = mod.code))
  
  
  
}