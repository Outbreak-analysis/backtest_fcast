library(R0)

source("idea.R")
source('../RESuDe_forecast/RESuDe_FCT.R')
source('../RESuDe_forecast/fit_stan.R')
source('../RESuDe_forecast/forecast.R')
source("../SEmInR/SEmInR_deterministic.R")
source("../SEmInR/SEmInR_deterministic_fit_mle.R")
source("../SEmInR/SEmInR_deterministic_fit_ABC.R")
source("../genGrowth/ggm_lib.R")

plot.GI <- function(g, GI.dist){
	n <- length(g)
	m <- sum(c(1:n)*g)
	v <- sum((c(1:n)^2)*g) - m^2
	plot(x=1:n, y=g, 
		 type="o",lwd=4,
		 pch=5,
		 xlab="time since infection",ylab="pdf",
		 main=paste0("GI distribution - ",GI.dist))
	abline(v=m, lwd=2, lty=2)
	abline(v=m+sqrt(v)*c(-1,1), lwd=1, lty=2)
	grid()
}

plot.exp.phase <- function(dat){
	# Check exponential phase is more or less present:
	loginc <- log(dat$inc+1)
	tt <-dat$t 
	m <- lm(loginc~tt)
	r2 <- summary(m)$r.squared
	plot(tt,loginc,
		 pch=16,type="o",
		 xlab="time",
		 main=paste0("Check exp phase: R2=",round(r2,3)))
	text(x = tt,y=loginc, labels = dat$inc, pos = 3)
	abline(a = m$coefficients[1], b=m$coefficients[2],lty=2)
	if( r2<0.5 || r2 > 0.91) 
		warning("Exponential phase either not satisfied, or too good (reaching peak incidence?)")
}

plot.fcast <- function(model,
					   dat,
					   inc.f.m,
					   inc.f.lo,
					   inc.f.hi,
					   Restim,
					   dat.full=NULL,
					   log=FALSE){
	
	tt <-dat$t 
	nf <- length(inc.f.m)
	ntt <- length(tt)
	f.rng <- (ntt+1):nf
	inc <- dat$inc
	ylab <- "incidence"
	
	if(log){
		inc.f.m <- log(inc.f.m)
		inc.f.lo <- log(inc.f.lo)
		inc.f.hi<- log(inc.f.hi)
		inc <- log(inc)
		ylab <- "LOG incidence"
	}
	title <- paste(model,"forecast time series\n",Restim)
	yrng <- c(0,max(inc.f.hi))
	if(!is.null(dat.full)) {
		inc.full <- dat.full$inc[f.rng]
		if(log) inc.full <- log(dat.full$inc[f.rng])
		yrng <- c(0,max(inc.f.hi,inc.full,na.rm = TRUE))
	}
	plot(x= c(tt,(length(tt)+1):nf),
		 y=inc.f.m,
		 cex=1, lwd=1, typ='o',
		 main = title,
		 xlab="time", ylab=ylab,
		 ylim=yrng)
	polygon(x=c(f.rng,rev(f.rng)),
			y = c(inc.f.lo[f.rng],rev(inc.f.hi[f.rng])),
			border = NA,
			col = rgb(0,0,0,0.2))
	points(inc,pch=16)
	
	if(!is.null(dat.full)){
		
		if(log) inc.full <- log(dat.full$inc[f.rng])
		
		points(x=dat.full$t[f.rng], 
			   y=inc.full,
			   pch=3, col="red",lwd=3)
		
		legend(x = "topleft",bty = "n",
			   pch=c(1,3), col=c("black","red"), 
			   pt.cex = 2,
			   legend = c("Forecast","Target"))
	}
	abline(v=tt[ntt],lty=2)
}

plot.fcast.vs.actual <- function(dat,
								 dat.full,
								 inc.f.m,
								 inc.f.lo,
								 inc.f.hi,
								 log=FALSE){
	tt <-dat$t 
	nf <- length(inc.f.m)
	ntt <- length(tt)
	f.rng <- (ntt+1):nf
	inc <- dat$inc
	inc.full <- dat.full$inc[f.rng]
	title <- "Incidence Forecast vs Actual"
	if(log){
		tiny <- 0.1 #10E-9
		inc.f.m <- log(inc.f.m+tiny)
		inc.f.lo <- log(inc.f.lo+tiny)
		inc.f.hi<- log(inc.f.hi+tiny)
		inc <- log(inc+tiny)
		inc.full <- log(inc.full+tiny)
		title <- "LOG Incidence Forecast vs Actual"
	}
	plot(x = inc.full,
		 y = inc.f.m[f.rng],
		 xlab = "Actual incidence", 
		 ylab = "Forecast incidence",
		 main = title,
		 ylim = range(inc.f.lo[f.rng],inc.f.hi[f.rng],inc.full,na.rm = TRUE),
		 cex=2,lwd=2)
	segments(x0=inc.full, x1=inc.full,
			 y0=inc.f.lo[f.rng], y1=inc.f.hi[f.rng],
			 lwd=1)
	grid()
	abline(0,1,col="red",lty=2,lwd=2)
}

plot.all <- function(model,
					 Restim,
					 dat,
					 dat.full,
					 inc.f.m,
					 inc.f.lo,
					 inc.f.hi){
	
	plot.fcast(model,
			   Restim,
			   dat = dat,
			   inc.f.m = inc.f.m,
			   inc.f.lo = inc.f.lo,
			   inc.f.hi = inc.f.hi,
			   dat.full = dat.full,
			   log=FALSE)
	plot.fcast.vs.actual(dat,
						 dat.full,
						 inc.f.m,
						 inc.f.lo,
						 inc.f.hi,
						 log=FALSE)
	plot.fcast(model,
			   Restim,
			   dat = dat,
			   inc.f.m = inc.f.m,
			   inc.f.lo = inc.f.lo,
			   inc.f.hi = inc.f.hi,
			   dat.full = dat.full,
			   log=TRUE)
	plot.fcast.vs.actual(dat,
						 dat.full,
						 inc.f.m,
						 inc.f.lo,
						 inc.f.hi,
						 log=TRUE)
	# Focus on the fit
	n2 <- nrow(dat)+1
	plot.fcast(model  = model,
			   Restim = Restim,
			   dat = dat,
			   inc.f.m  = inc.f.m[1:n2],
			   inc.f.lo = inc.f.lo[1:n2],
			   inc.f.hi = inc.f.hi[1:n2],
			   dat.full = NULL,
			   log=FALSE)
}


translate.model <- function(x){
	res <- NA
	if(x=="WalLip") res <- "EG"
	if(x=="WhiPag") res <- "ML"
	if(x=="SeqBay") res <- "SB"
	if(x=="CoriParam") res <- "ParametricSI"
	if(x=="CoriNonParam") res <- "NonParametricSI"
	if(x=="CoriUncertain") res <- "UncertainSI"
	if(x=="IDEA") res <- "IDEA"
	
	if(is.na(res)) stop(paste("Model name unknown:",x))
	return(res)
}

unpack.prm <- function(prms){
	### Unpack parameters depending on model:
	### WARNING: in GLOBAL environment
	
	model <<- prms[["model"]]
	
	pname <- c("dat","dat.full", "horiz.fcast",
			   "GI.dist","GI.val","GI.truncate")
	
	pname.cori <- c("cori.window","cori.mean.prior","cori.std.prior",
					"Std.Mean.SI","Min.Mean.SI",
					"Max.Mean.SI","Std.SI", "Std.Std.SI", 
					"Min.Std.SI", "Max.Std.SI", 
					"n.coriUnc.meanstdv", "n.coriUnc.postR")
	
	pname.resude <- c('pop_size', 'GI_span', 'GI_mean', 'GI_var',
					  'alpha','kappa',
					  'mcmc_iter', 'mcmc_nchains','mcmc_diagnostic')
	
	pname.seminr <- c('prm.fxd','prm.to.fit')
	
	if(grepl("Cori",model)) pname <- c(pname, pname.cori)
	if(model=='RESuDe')     pname <- c(pname, pname.resude)
	if(model=='SEmInRdet')  pname <- c(pname, pname.seminr)
	
	for(p in pname) assign(x = p,
						   value = prms[[p]],
						   envir = .GlobalEnv)
}


### - - - - - - - - - - - - - - - - -
### - - - - FITTING FUNCTIONS - - - - 
### - - - - - - - - - - - - - - - - -

fit.renewal <- function(prms){
	### FIT THE REPRODUCTIVE NUMBER
	### USING 'R0' and'EpiEstim' PACKAGES
	
	unpack.prm(prms)
	
	model2 <- translate.model(model)
	
	### SWITCH BETWEEN REQUESTED MODELS 
	
	GI <- NULL
	
	# R0 package:
	if(model %in% c("WalLip","WhiPag","SeqBay")){
		# Create generation time 
		GI <- generation.time(type     = GI.dist, 
							  val      = GI.val,
							  truncate = GI.truncate,
							  step     = 1 )
		# Estimate R
		inc <- dat$inc
		if (model=="SeqBay"){
			# There is a conceptual problem with this model
			# when incidence data = 0, because:
			# I[t+1] ~ Poisson(*I[t])
			# so when I[t]=0, probability(I[t+1]>0)=0
			# (bc the poisson intensity is 0)
			#
			# So, try to go round this issue by
			# artificially nudging incidence away from 0 (e.g. to 1)
			if(any(inc==0)) {
				warning("Incidence was changed for SeqBay method")
				inc[inc==0] <- 1
			}
		}
		R <- estimate.R(epid    = inc,
						GT      = GI, 
						methods = model2)
	}
	
	# EpiEstim package
	if(grepl(pattern = "Cori",x = model)){
		### Reference: Cori 2013 American Journal of Epidemiology
		# We just want to estimate R at the 
		# last incidence point available:
		# (this model provides a 'sliding' estimation for R)
		N <- length(dat$inc)
		# deal with default values:
		if(is.null(cori.mean.prior)) cori.mean.prior <- 5
		if(is.null(cori.std.prior)) cori.std.prior <- 5
		
		R <- EstimateR(I       = dat$inc,
					   T.Start = (N-cori.window+1),  # <-- Estimate R at the end of the data set (bc forecast)
					   T.End   = N,                  # <-- Estimate R at the end of the data set (bc forecast)
					   method  = model2,
					   #SI.Distr = NULL,
					   Mean.SI    = GI.val[1],
					   Std.SI     = GI.val[2],
					   Mean.Prior = cori.mean.prior,
					   Std.Prior  = cori.std.prior,
					   plot = FALSE,
					   # Param below for "CoriUncertain" model only:
					   Std.Mean.SI = Std.Mean.SI,
					   Min.Mean.SI = Min.Mean.SI, 
					   Max.Mean.SI = Max.Mean.SI, 
					   Std.Std.SI  = Std.Std.SI, 
					   Min.Std.SI  = Min.Std.SI, 
					   Max.Std.SI  = Max.Std.SI, 
					   n1 = n.coriUnc.meanstdv, 
					   n2 = n.coriUnc.postR
		)
	}
	### Retrieve estimates depending on model chosen:
	if(model %in% c("WalLip","WhiPag")){
		R.m  <- R$estimates[[model2]]$R
		R.lo <- R$estimates[[model2]]$conf.int[1]
		R.hi <- R$estimates[[model2]]$conf.int[2]
	}
	if(model %in% c("SeqBay")){
		R.m <- R$estimates[[model2]]$R
		# Take last estimate:
		nr   <- length(R.m)
		R.m  <- R.m[nr]
		R.lo <- R$estimates[[model2]]$conf.int[nr,1]
		R.hi <- R$estimates[[model2]]$conf.int[nr,2]
	}
	if(grepl("Cori",model)){
		R.m  <- R$R$`Mean(R)`
		R.lo <- R$R$`Quantile.0.025(R)`
		R.hi <- R$R$`Quantile.0.975(R)`
		R.SIDistr <- R$SIDistr
	}
	return(list(R.m = R.m, 
				R.lo = R.lo, 
				R.hi = R.hi,
				R = R,
				GI = GI,
				model = model))
}

### Fit generailzed growth model (Chowell et al 2016)
fit.GGM <- function(prms){
	
	inc <- dat$inc
	n   <- length(inc)
	dat <- data.frame(t=1:length(inc), inc=inc)
	# guess for fitted param:
	prm.init <- c(r = (inc[n]-inc[1])/n, 
				  p = 0.9)
	
	fit <- estimate.CI(dat          = dat, 
					   CIwidth      = 0.95,
					   n.MC         = 200,
					   prm.init     = prm.init, 
					   relative.err = FALSE)
	return(fit)
}


fit.seqBay <- function(prms){
	### FIT REPRODUCTIVE NUMBER 
	### FOR SEQUENTIAL BAYESIAN METHOD
	###
	unpack.prm(prms)
	stopifnot(model=='SeqBay')
	model2 <- translate.model(model)
	GI <- NULL
	inc <- dat$inc
	# R0 package:
	# Create generation time 
	GI <- generation.time(type     = GI.dist, 
						  val      = GI.val,
						  truncate = GI.truncate,
						  step     = 1 )
	# Estimate R:
	# There is a conceptual problem with this model
	# when incidence data = 0, because:
	# I[t+1] ~ Poisson(*I[t])
	# so when I[t]=0, probability(I[t+1]>0)=0
	# (bc the poisson intensity is 0)
	#
	# So, try to go round this issue by
	# artificially nudging incidence away from 0 (e.g. to 1)
	if(any(inc==0)) {
		warning("Incidence was changed for SeqBay method")
		inc[inc==0] <- 1
	}
	
	R <- estimate.R(epid    = inc,
					GT      = GI, 
					methods = model2)
	
	R.m <- R$estimates[[model2]]$R
	# Take last estimate:
	nr   <- length(R.m)
	R.m  <- R.m[nr]
	R.lo <- R$estimates[[model2]]$conf.int[nr,1]
	R.hi <- R$estimates[[model2]]$conf.int[nr,2]
	
	return(list(R.m  = R.m, 
				R.lo = R.lo, 
				R.hi = R.hi,
				R    = R,
				GI   = GI,
				model = model))
}

fit.resude <- function(prms) {
	
	# unpack paramters:
	unpack.prm(prms)
	
	# forecasting horizon:
	dat.obs  <- dat$inc
	last.obs <- length(dat.obs)
	
	# effective population bounds:
	pop_mean <-  resude.pop.init
	pop_lsd  <-  2.0
	pop_hi <- qlnorm(p=0.99, meanlog = log(pop_mean),sdlog = pop_lsd)
	pop_lo <- qlnorm(p=0.1, meanlog = log(pop_mean),sdlog = pop_lsd)
	pop_lo <- max(pop_lo,sum(dat.obs)+1)# <-- pop cannot be smaller than cumul incidence!
	pop_hi <- round(pop_hi,0)
	pop_lo <- round(pop_lo,0)
	
	# Define Stan's known data:
	data.stan <- list(numobs   = last.obs,
					  Iobs     = dat.obs,
					  R0_lo    = 0.7,
					  R0_hi    = 10,
					  GI_span  = GI_span,
					  GI_mean  = GI_mean,
					  GI_var   = GI_var,
					  alpha    = 0,
					  kappa    = 0,
					  pop_hi   = pop_hi,
					  pop_lo   = pop_lo, 
					  pop_mean = pop_mean,
					  pop_lsd  = pop_lsd # log stddev
	)
	print("RESuDe parameter before Stan fit:")
	print(data.stan)
	
	# Fit RESuDe model using Stan:
	FIT <- RESuDe.fit.stan(model.filename = 'fit-resude-light.stan', 
						   dat      = data.stan, 
						   n.iter   = mcmc_iter, 
						   n.chains = mcmc_nchains,
						   #n.cores = floor(parallel::detectCores()/2),
						   plot.compTruth = FALSE
	) 
	# Show diagnostic plots for Stan fit:
	if(mcmc_diagnostic){
		np <- names(FIT$prm.sample)
		np <- np[np!="Iout"]
		pairs(FIT$fit)	
		np <- np[np!="lp__"]
		traceplot(FIT$fit, pars=np, alpha=0.5,inc_warmup=TRUE)
	}
	return(FIT)
}

rough.R0 <- function(inc.obs){
	# Uses R0 ~ 1 / sum(exp(-r*t)*g(t)) 
	n <- length(inc.obs)
	t.obs <- 1:n # assumes regular observations
	loginc <- log(inc.obs)
	lminc <- lm(loginc~t.obs)
	r <- lminc$coefficients[2]
	ngi <- min(n,1)
	gi <- 1/ngi # uniform generation interval distribution :-\
	roughR0 <- sum(exp(-r*t.obs[1:ngi])*gi)
	return(roughR0)
}


fit.SEmInRdet <- function(prms, fit.type){
	
	# unpack paramters:
	unpack.prm(prms)
	inc.obs  <- dat$inc
	t.obs    <- 1:length(inc.obs)
	
	
	if(fit.type == 'mle'){
		print('DEBUG: SEmInR mle fit with inital values:')
		print(prm.to.fit)
		logparam <- TRUE
		if(logparam) prm.to.fit <- log(prm.to.fit)
		# minimization:
		FIT <- fit.mle.SEmInR(prm.to.fit, 
							  prm.fxd, 
							  t.obs, 
							  inc.obs,
							  logparam = logparam,
							  method = "Nelder-Mead", #"SANN",# "CG",# "Nelder-Mead",#'SANN',#"L-BFGS-B",
							  maxit = 500)
		
		prm.fitted <- FIT[['prm.fitted']]
		llkmin     <- FIT[['llkmin']]    
		print('DEBUG: SEmInR mle fitted values:')
		print(prm.fitted)
		
		# Try to approximate CI value:
		cival <- CI.llk.sample(CIlevel = 0.95, 
							   nsample = 100, 
							   prm.fitted, 
							   llkmin, prm.fxd, t.obs,
							   inc.obs = inc.obs,
							   prop.search = 0.5)
		M <- matrix(unlist(cival), ncol=length(prm.fitted),byrow = T)
		if("R0" %in% names(prm.fitted)){
			j <- which(names(prm.fitted)=="R0")
			R0.lo <- min(M[,j])
			R0.hi <- max(M[,j])
			R0    <- prm.fitted['R0']
		}
		return(list(fit = FIT, cival = cival, 
					R0=R0, R0.lo=R0.lo, R0.hi=R0.hi))
	}
	if(fit.type == 'ABC'){
		### TO DO : READ FROM A FILE !!!
		priors <- list(infectious_mean = list('unif', prm=list(1,9)),
					   latent_mean     = list('unif', prm=list(1,9)),
					   popSize         = list('unif', prm=list(1e3,1e6)),
					   R0              = list('unif', prm=list(0.95,9))
		)
		priors.prm.to.fit <- priors
		
		FIT <- fit.ABC.SEmInR(t.obs, inc.obs, 
							  prm.fxd, priors.prm.to.fit, 
							  nABC = 1E4, post.prop=0.01)
		return(FIT)
	}
}


### - - - - - - - - - - - - - - - - -
### - - - - SIMULATION FUNCTIONS - - - - 
### - - - - - - - - - - - - - - - - -

setup_GI_renewal <- function(model,GI, R){
	g <- NULL
	if(model %in% c("WalLip","WhiPag","SeqBay")){
		#  Remove day 0 of GI distribution
		if(GI$time[1]==0) g <- GI$GT[-1]  
	}
	if(grepl("Cori",model)){
		if(model=="CoriNonParam") {
			g <- GI.val
			GI.dist <- "empirical"
		}
		if(model=="CoriParam") {
			g <- R$SIDistr
			GI.dist <- "gamma"
		}
		if(model=="CoriUncertain") 
			stop("Model CoriUncertain implementation not finished!")
		
		#  Remove day 0 of GI distribution
		if(g[1,1]==0) g <- g[-1,] 
		g <- g[,2]
	}
	return(list(g=g, GI.dist=GI.dist))
}

simFwd_renewal <- function(obsinc, # observed incidence
						   g,
						   R.m,R.lo,R.hi,
						   horiz.fcast){
	
	inc.f.m <- inc.f.lo <- inc.f.hi <- obsinc
	n.i <- length(inc.f.m)
	
	# Apply renewal equation at each time step 
	# until forecast horizon:
	for(i in 1:horiz.fcast){
		n.g2 <- min(length(g), n.i)
		inc.rev.m <- rev(inc.f.m)[1:n.g2]
		inc.rev.lo <- rev(inc.f.lo)[1:n.g2]
		inc.rev.hi <- rev(inc.f.hi)[1:n.g2]
		
		# renewal equation:
		inc.f.m  <- c(inc.f.m,  R.m  * sum(g[1:n.g2] * inc.rev.m[1:n.g2]))
		inc.f.lo <- c(inc.f.lo, R.lo * sum(g[1:n.g2] * inc.rev.lo[1:n.g2]))
		inc.f.hi <- c(inc.f.hi, R.hi * sum(g[1:n.g2] * inc.rev.hi[1:n.g2]))
	}
	return(list(inc.f.m=inc.f.m,
				inc.f.lo=inc.f.lo,
				inc.f.hi=inc.f.hi))
}

simFwd_seqBay <- function(obsinc, # observed incidence
						  GI,
						  R.m,R.lo,R.hi,
						  horiz.fcast){
	# Apply the formula I(t+h) = exp(h*gamma(R-1))*I(t)
	# from Bettencourt 2008 PLoS ONE
	n.i <- length(obsinc)
	tmp.m  <- exp( (1:horiz.fcast)*(R.m-1.0)/GI$mean )
	tmp.lo <- exp( (1:horiz.fcast)*(R.lo-1.0)/GI$mean )
	tmp.hi <- exp( (1:horiz.fcast)*(R.hi-1.0)/GI$mean )
	inc.f.m  <- c(obsinc, obsinc[n.i] * tmp.m)
	inc.f.lo <- c(obsinc, obsinc[n.i] * tmp.lo)
	inc.f.hi <- c(obsinc, obsinc[n.i] * tmp.hi)
	
	return(list(inc.f.m  = inc.f.m,
				inc.f.lo = inc.f.lo,
				inc.f.hi = inc.f.hi))
}

simFwd_GGM <- function(fit, horiz.fcast,obsinc){
	
	tvec <- 1:(length(obsinc)+horiz.fcast)
	sim.md <- genGrowth.inc(c0=obsinc[1],
							r = fit[['r.md']],
							p = fit[['p.md']],
							tvec = tvec)
	
	sim.lo <- genGrowth.inc(c0=obsinc[1],
							r = fit[['r.ci']][1],
							p = fit[['p.ci']][1],
							tvec = tvec)
	
	sim.hi <- genGrowth.inc(c0=obsinc[1],
							r = fit[['r.ci']][2],
							p = fit[['p.ci']][2],
							tvec = tvec)
	
	return(list(inc.f.m  = sim.md,
				inc.f.lo = sim.lo,
				inc.f.hi = sim.hi))
}

simFwd_RESuDe <- function(FIT, 
						  horiz.fcast,
						  GI_span, GI_mean, GI_var, 
						  alpha,
						  kappa){
	
	obs.inc <- FIT$prm.sample$Iout[1,]
	last.obs <- length(obs.inc)
	
	FCAST <- RESuDe.forecast.light(prm = FIT$prm.sample, # <-- sampled parameter values after the fit
								   fcast.horizon = horizon,
								   alpha = alpha,
								   kappa = kappa, 
								   GI_span = GI_span, 
								   GI_mean = GI_mean,
								   GI_var  = GI_var,
								   last.obs = last.obs,
								   nmax = 500,
								   syn.inc.full = NULL, 
								   do.plot = FALSE,
								   CI1 = 50, CI2 = 95,
								   seed=123)
	
	inc.f.m  <- FCAST$fcast.cone[,3]
	inc.f.lo <- FCAST$fcast.cone[,1]
	inc.f.hi <- FCAST$fcast.cone[,5]
	
	nobs <- length(obs.inc)
	nfor <- length(inc.f.m)
	inc.f.m  <- c(obs.inc, inc.f.m[ (nobs+1):nfor])
	inc.f.lo <- c(obs.inc, inc.f.lo[(nobs+1):nfor])
	inc.f.hi <- c(obs.inc, inc.f.hi[(nobs+1):nfor])
	
	return(list(inc.f.m  = inc.f.m,
				inc.f.lo = inc.f.lo,
				inc.f.hi = inc.f.hi,
				inc.f.all = FCAST$sf))
}


simFwd_SEmInRdet <- function(prm.fitted, prm.fxd, cival) {
	
	# Simulate forward with fitted data (point estimate):
	sim.fit  <- simul.SEmInR.det(prm.fitted, prm.fxd)
	df.fit   <- sim.fit$ts
	inc.best0 <- df.fit$inc
	
	t <- df.fit$time
	dt <- df.fit$time[2]-df.fit$time[1]
	tmax <- max(t)
	idx <- which(abs(t-round(t))<dt/2)
	if(t[1]==0) idx <- idx[-1] # remove 0
	
	# Simulate forward with fitted data (CI envelop):
	inc.CI <- list()
	for(s in 1:length(cival)){
		if(s%%10==0) print(paste('simulating CI incidence',s,'/',length(cival)))
		sim.CI <- simul.SEmInR.det(cival[[s]],prm.fxd)
		inc.CI[[s]] <- sim.CI$ts$inc
	}
	
	m.CI <- matrix(unlist(inc.CI),ncol=length(inc.CI[[1]]),byrow = TRUE)
	inc.ci.lo0 <- apply(m.CI,MARGIN = 2,FUN=min)
	inc.ci.hi0 <- apply(m.CI,MARGIN = 2,FUN=max)
	
	# Take values at integer times:
	inc.best  <- inc.best0[idx]
	inc.ci.lo <- inc.ci.lo0[idx]
	inc.ci.hi <- inc.ci.hi0[idx]
	
	return(list(inc.f.m  = inc.best,
				inc.f.lo = inc.ci.lo,
				inc.f.hi = inc.ci.hi
	))
}

### - - - - - - - - - - - - - - - - -
### - - - -   FORECASTS   - - - - 
### - - - - - - - - - - - - - - - - -

### FORECAST INCIDENCE WITH A CHOICE OF MODELS
###
fcast_incidence <- function(prms, do.plot=FALSE){
	
	unpack.prm(prms)
	
	print(paste("DEBUG::",model))
	
	### - - - - - - - - - - - - 
	###   Perform fit to data
	### - - - - - - - - - - - - 
	
	if(model %in% c("WalLip","WhiPag") || grepl("Cori",model)){
		fit <- fit.renewal(prms)	
	} 
	if(model %in% "SeqBay") fit <- fit.seqBay(prms)
	if(model == 'GGM')      fit <- fit.GGM(prms)
	if(model == "RESuDe")   fit <- fit.resude(prms)
	if(model == "SEmInRdet") {
		seminr.fit.type <- 'mle'  # 'mle' or 'ABC'
		fit <- fit.SEmInRdet(prms,fit.type = seminr.fit.type)
	}
	
	### Retrieve fit results:
	
	if(model %in% c("SeqBay","WalLip","WhiPag")|| grepl("Cori",model)){
		R.m   <- fit[["R.m"]]
		R.lo  <- fit[["R.lo"]]
		R.hi  <- fit[["R.hi"]]
		R     <- fit[["R"]]
		GI    <- fit[["GI"]]
		
		gitmp   <- setup_GI_renewal(model, GI, R)
		GI.dist <- gitmp[["GI.dist"]]
		g       <- gitmp[["g"]]
	}
	if(model == "RESuDe"){
		R    <- fit$prm.sample$R0
		R.lo <- quantile(R,probs = 0.5-CI/2)
		R.hi <- quantile(R,probs = 0.5+CI/2)
		R.m  <- median(R)
	}
	if(model == 'GGM'){
		# may have to change that; place holder for now...
		R <- R.lo <- R.hi <- R.m <- fit[['r.md']]
	}
	if(model == 'SEmInRdet'){
		if(seminr.fit.type=='mle'){
			R.lo <- fit[['R0.lo']]
			R.hi <- fit[['R0.hi']]
			R.m  <- fit[['R0']]
			R    <- R.m	
		}
		if(seminr.fit.type=='ABC'){
			R0post <- fit[['posteriors']][,'R0']
			R.lo <- quantile(R0post,probs = 0.5-CI/2)
			R.m  <- quantile(R0post,probs = 0.5)
			R.hi <- quantile(R0post,probs = 0.5+CI/2)
			R    <- R.m
		}
	}
	
	### Simulate forward (with fitted model parameters)
	
	inc.f.m <- inc.f.lo <- inc.f.hi <- dat$inc
	n.i <- length(inc.f.m)
	
	if(model %in% c("WalLip","WhiPag") || grepl("Cori",model)){
		sim <- simFwd_renewal(obsinc = dat$inc, # observed incidence
							  g,
							  R.m, R.lo, R.hi,
							  horiz.fcast)
	}
	if(model %in% c("SeqBay")){
		sim <- simFwd_seqBay(obsinc = dat$inc, 
							 GI, R.m, R.lo, R.hi, horiz.fcast)
	}
	
	if (model=="IDEA") {
		### IDEA impementation fit 
		### and forecast at the same time
		sim <- idea.forecast(data = dat$inc,
							 stoch = F,
							 CI = 0.999,
							 horiz.fcast = horiz.fcast,
							 GI = GI.val[1],
							 ignore.data = 1)
		
		R <- R.m <- R.hi <- R.lo <- sim[["R0"]]
	}
	if (model == 'GGM'){
		sim <- simFwd_GGM(fit, horiz.fcast, obsinc = dat$inc)
	}
	
	if(model=='RESuDe'){
		sim <- simFwd_RESuDe(FIT = fit,
							 horiz.fcast = horiz.fcast,
							 GI_span = GI_span,
							 GI_mean = GI_mean, 
							 GI_var  = GI_var,
							 alpha   = alpha, 
							 kappa   = kappa)
	}
	
	if(model == "SEmInRdet"){
		if(seminr.fit.type=='mle'){
			prm.fitted <- fit$fit[['prm.fitted']]
			sim <- simFwd_SEmInRdet(prm.fitted,prm.fxd,fit$cival)
		}
		if(seminr.fit.type=='ABC'){
			postinc <- post.incidence(fit[['posteriors']], CI)
			# simulations returned are with a different time
			# granularity (bc ODE solved), so extract only at rounded times:
			horiz <- prms$prm.fxd['horizon']
			sim.time <- postinc$time
			sim <- list()
			sim[["inc.f.m"]]   <- time.subset(t=sim.time, x=postinc[['inc.md']],t.subset = 1:horiz)
			sim[["inc.f.lo"]]  <- time.subset(t=sim.time, x=postinc[['inc.lo']],t.subset = 1:horiz)
			sim[["inc.f.hi"]]  <- time.subset(t=sim.time, x=postinc[['inc.hi']],t.subset = 1:horiz)
			
			# sim[["inc.f.m"]]  <- postinc[['inc.md']]
			# sim[["inc.f.lo"]] <- postinc[['inc.lo']]
			# sim[["inc.f.hi"]] <- postinc[['inc.hi']]
		}
	}
	
	### Retrieve all simulated incidences:
	inc.f.m  <- sim[["inc.f.m"]]
	inc.f.lo <- sim[["inc.f.lo"]]
	inc.f.hi <- sim[["inc.f.hi"]]
	
	# R estimates:
	Restim <- paste0(round(R.m,2)," (",
					 round(R.lo,2),";",
					 round(R.hi,2),")")
	
	#  - - - Plots - - -
	if (do.plot){
		par(mfrow=c(2,3))
		try(plot.exp.phase(dat),silent = T)
		# try(plot.GI(g,GI.dist),silent = T)
		try(plot.all(model, Restim,dat, dat.full, 
					 inc.f.m, inc.f.lo, inc.f.hi),
			silent = F)
	}
	# Target data (if exists, for testing purpose)
	target.dat <- NULL
	if(!is.null(dat.full)) target.dat <- dat.full$inc[(length(dat$t)+1):length(inc.f.m)]
	
	return(list(R=R,
				inc.f.m = inc.f.m,
				inc.f.lo = inc.f.lo,
				inc.f.hi = inc.f.hi,
				obs.dat = dat$inc,
				target.dat = target.dat) )
}

