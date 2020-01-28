model
{
  #### counts and overdispersion effects  ###### Link-Sauer model with nested observer-route effect
  for( k in 1 : ncounts )
  {
    Elambdas[k] <- beta[strat[k]] * (year[k] - fixedyear) + eta*firstyr[k] + gam.sm[k] + obs[strat[k],obser[k]] + strata[strat[k]] + yeareffect[year[k],strat[k]]
    Elambda[k] ~ dnorm(Elambdas[k], taunoise)
    log(lambda[k]) <- Elambda[k]
    
    count[k] ~ dpois(lambda[k])
    #----------------------------------#
    fcount[k] ~ dpois(lambda[k])
    err[k] <- pow(count[k]-lambda[k],2)/lambda[k]
    ferr[k] <- pow(fcount[k]-lambda[k],2)/lambda[k]
    fzero[k] <- equals(fcount[k],0)
    loglik[k] <- logdensity.pois(count[k], lambda[k])
    #----------------------------------#
  }
  
  nfzero <- sum(fzero[1:ncounts])
  gof <- sum(err[1:ncounts])
  fgof <- sum(ferr[1:ncounts])
  diffgof <- gof-fgof
  posdiff <- step(diffgof)
  
  taunoise ~ dgamma(0.001,0.001)
  sdnoise <- 1 / pow(taunoise, 0.5)
  sdobs <- 1 / pow(tauobs, 0.5)
  tauobs ~ dgamma(0.001,0.001)
  eta ~ dnorm( 0.0,1.0E-6)
  STRATA ~ dnorm( 0.0,0.01)
  taustrata ~ dgamma(0.001,0.0001) #<- 1/pow(sdbeta,2)#
  sdstrata <- 1/pow(taubeta,0.5)#~ dunif(0.001,10)
  
  #---------------------------------------------------------#
  
  #----------------------------------#
  #### stratum effects  ######
  for( s in 1 : nstrata )
  {
    #### observer effects  ######
    
    for( i in 1 : nobservers[s] )
    {
      obs[s,i] ~ dnorm( 0.0,tauobs)
    }
    
    #### end observer effects  ######
    
    beta[s] ~ dnorm(BETA,taubeta)
    strata[s] ~ dnorm(STRATA,taustrata)
    sdyear[s] <- 1/pow(tauyear[s],0.5)
    tauyear[s] ~ dgamma(0.001,0.001) #
    
    #### stratum specific year effects  ######
    for( y in ymin : ymax )
    {
      yeareffect[y,s] ~ dnorm( 0.0, tauyear[s])
    }
    
    expstrata[s] <- exp(strata[s])
    overdisp[s] <- 1 + 1/(expstrata[s]*taunoise)
    #-------------------------------------------------#
    
  }# end s strata loop
  
  ################################## seasonal GAM
  
  for ( i in 1:ncounts ){
    for(k in 1:nknots){
      X.part[i,k] <- beta.X[strat[i],k]*(X.basis[i,k])   #######GAM smooth component
    }
    gam.sm[i] <- sum(X.part[i,1:nknots])
  }#i
  
  
  ###################################################################
  ##### this text assumes the GAM parameters (betas) are estimated as random effects########
  ##### your will need to provide the variable ngroups (an integer indicating the number of groups in the random effect)####
  ##### and ensure that your R-script retains the correct order of groups in relation to the predicted values indexed by the 1:ngroups section below ######
  
  # Following Crainiceanu, C. M., Ruppert, D. & Wand, M. P. (2005). Bayesian Analysis for Penalized Spline Regression Using WinBUGS. Journal of Statistical Softare, 14 (14), 1-24.
  
  tauX~dgamma(1.0E-2,1.0E-4) #alternate prior, original from Cainiceanu et al. second gamma parameter == 0.0001 << (abs(mean(B.X[]))^2)/2, mean(B.X[]) ~ 0.2
  sdX <- 1/(pow(tauX,0.5)) # 
  taubetagam <- 1/pow(sdbetagam,2) # prior on precision of gam coefficients(
  sdbetagam ~ dunif(0,5)
  
  for(k in 1:nknots)
  { # Computation of GAM components
    B.X[k] ~ dnorm(0,tauX)
    
    for(j in 1:nstrata)
    {
      beta.X[j,k] ~ dnorm(B.X[k],taubetagam) ####B hyperparameter for that knot k
      
      
      for ( i in 1:npredpoints ){ #this allows the user to model just the GAM component at particular values of the predictor
        X.partpred[i,k,j] <- beta.X[j,k]*(X.basispred[i,k])
        
      }#i
      
    }#k
  }#j
  
  for (i in 1:npredpoints){    #predpoints are dates in the season
    for(j in 1:nstrata)
    {
      gam.smooth[i,j] <- sum(X.partpred[i,1:nknots,j])  ###### take a look at the different smooth function per strata, what is the seasonal effect
    }#k
  }#i
  
  
  ####################### end seasonal GAM
  
  #### hierarchical trend effects  ######
  
  
  BETA ~ dnorm(0,0.01)
  
  taubeta ~ dgamma(1,0.0001) #
  ### informative taubeta, parameterized so that:
  ### mean of taubeta = ~10000
  ### the mean of sdbeta =~ 0.02
  ### there is <5% prob of sdbeta > 0.04 (i.e., there is a small, but non-zero probability that the SD of the stratum beta
  ### terms is large enough to to encompass trends that vary from -8%/year to 8%/year with a continental trend = 0)
  ### this informative prior should only be useful when there are few strata
  
  sdbeta <- 1/pow(taubeta,0.5)#
  
  #for( i in 1 : nstrata )
  #{
  #	for( t in ymin : ymax )
  #	{
  #		n[i,t] <- nonzeroweight[i]*exp(strata[i] + beta[i]*(t-fixedyear) + yeareffect[t,i] + 0.5*sdnoise*sdnoise+ 0.5*sdobs*sdobs)
  #	}
  #}
  
  #### summary statistics  ######
  sdn <- exp(0.5*sdnoise*sdnoise)
  
  for( i in 1 : nstrata )
  {
    for( t in ymin : ymax )
    {
      for(o in 1 : nobservers[i])
      {
        no[i,t,o] <- exp(strata[i] + beta[i]*(t-fixedyear) + yeareffect[t,i] + obs[i,o])
      }
      
      mn[i,t] <- mean(no[i,t,1 : nobservers[i]])
      n[i,t] <- nonzeroweight[i]*(mn[i,t]*sdn)
    }
  }
  
  #-------------------------------------------------#
}

