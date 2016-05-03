
model.colour <- function(model.name){
	col <- "grey"
	if(model.name == "Cori") col <- "black"
	if(model.name == "WalLip") col <- "blue"
	if(model.name == "WhiPag") col <- "green"
	if(model.name == "SeqBay") col <- "red"
	if(model.name == "IDEA") col <- "orange"
	return(col)
}

compare.fcast.early.shiny <- function(fcast, dolog){
	### Compare forecasts for Shiny GUI
	
	n.model <- length(fcast)
	tgt <- fcast[[1]]$target.dat
	n.tgt <- length(tgt)
	n <- length(fcast[[1]]$inc.f.m)
	frng <- (n-n.tgt+1):n
	n.obs <- length(fcast[[1]]$obs.dat)
	
	obs.dat <- fcast[[1]]$obs.dat
	tgt.dat <- fcast[[1]]$target.dat
	
	ymax <- 0
	for (i in 1:length(fcast)) 
		ymax <- max(ymax, fcast[[i]]$inc.f.hi, tgt.dat)
	
	ymin <- ifelse(dolog,1,0)
	plot(0, cex=0,
		 xlim = c(1,length(fcast[[i]]$inc.f.m)),
		 ylim = c(ymin,ymax),
		 las = 1,
		 log = ifelse(dolog,"y",""),
		 xlab = "time",
		 ylab = "")
	grid(lty=2,col="lightgrey")
	
	nudge <- 0.15*c(0:(length(fcast)-1))
	nudge <- nudge-mean(nudge)
	
	# Plot observed data:
	points(x = 1:n.obs, y = obs.dat, typ="s")
	points(x = 1:n.obs, y = obs.dat, typ="p", pch=16)
	# Plot future data if available:
	if(!is.null(tgt.dat)){
		points(x = frng, y = tgt.dat, 
			   cex = 3,
			   col = "lightgrey",
			   pch = 15)
	}
	lwd.fcast <- 4
	
	for (i in 1:n.model) {
		f.m <- fcast[[i]]$inc.f.m[frng]
		f.hi <- fcast[[i]]$inc.f.hi[frng]
		f.lo <- fcast[[i]]$inc.f.lo[frng]
		
		points(x = frng+nudge[i],
			   y = f.m,
			   col = model.colour(names(fcast[i])), 
			   lwd = lwd.fcast, 
			   cex = 1.3)
		segments(x0 = frng+nudge[i], 
				 x1=frng+nudge[i],
				 y0 = f.lo, y1 = f.hi, 
				 col = model.colour(names(fcast[i])), 
				 lwd = lwd.fcast/2)
	}
	legend(x="topleft",
		   legend = names(fcast),
		   col = sapply(names(fcast),FUN = model.colour),
		   lwd = lwd.fcast, 
		   pch = 1)
}

create.model.prm <- function(dat,
							 dat.full,
							 horizon,
							 horiz.fcast ,  
							 GI.mean,GI.stdv,
							 GI.dist,
							 cori.window,
							 GI_span,
							 pop_size,
							 mcmc_iter,
							 mcmc_nchains,
							 mcmc_diagnostic,
							 SEmInR.prm.fxd,
							 SEmInR.prm.to.fit){
	
	### CREATES A LIST WITH ALL MODELS USED FOR FORECASTING
	### Add/remove models here
	
	PRM <- list(Cori = list(model = "CoriParam",  
							dat = dat,
							dat.full = dat.full,
							horiz.fcast = horiz.fcast,  
							GI.val = c(GI.mean,GI.stdv),
							cori.window = cori.window)
				,
				WalLip = list(model = "WalLip",
							  dat = dat,
							  dat.full = dat.full,
							  horiz.fcast = horiz.fcast,  
							  GI.dist = "gamma",  # gamma, lognormal, weibull
							  GI.val = c(GI.mean,GI.stdv))
				,
				WhiPag = list(model = "WhiPag",
							  dat = dat,
							  dat.full = dat.full,
							  horiz.fcast = horiz.fcast,  
							  GI.dist = "gamma",  # gamma, lognormal, weibull
							  GI.val = c(GI.mean,GI.stdv))
				,
				SeqBay = list(model = "SeqBay",
							  dat = dat,
							  dat.full = dat.full,
							  horiz.fcast = horiz.fcast,  
							  GI.dist = "gamma",  # gamma, lognormal, weibull
							  GI.val = c(GI.mean,GI.stdv))
				,
				IDEA = list(model = "IDEA",
							dat = dat,
							dat.full = dat.full,
							horiz.fcast = horiz.fcast,  
							GI.val = GI.mean)
				,
				RESuDe = list(model = "RESuDe",
							  dat = dat,
							  dat.full = dat.full,
							  horiz.fcast = horiz.fcast,
							  horizon = horizon,
							  GI_span = GI_span,
							  pop_size = pop_size,
							  mcmc_iter = mcmc_iter,
							  mcmc_nchains = mcmc_nchains,
							  mcmc_diagnostic = mcmc_diagnostic)
				,
				SEmInRdet = list(model = "SEmInRdet",
								 dat = dat,
								 dat.full = dat.full,
								 horiz.fcast = horiz.fcast,
								 prm.fxd = SEmInR.prm.fxd,
								 prm.to.fit = SEmInR.prm.to.fit
				)
	)
	return(PRM)
}

fcast.wrap.mc <- function(m, dat.all,
						  ttrunc,
						  horizon,
						  horiz.fcast,
						  GI.mean,GI.stdv,
						  cori.window,
						  GI_span,
						  pop_size,
						  mcmc_iter,
						  mcmc_nchains,
						  mcmc_diagnostic,
						  SEmInR.prm.fxd,
						  SEmInR.prm.to.fit,
						  rel.err,
						  do.plot,
						  single.model.fcast = NULL){
	
	### WRAP FOR PARALLEL CODE IN THE MONTE CARLO LOOP
	
	# Extract one realization:
	dat.full <- subset(dat.all, mc == m)
	# Find practical start of epidemic:
	tstart <- find.epi.start(dat.full$inc, w=4, thres.rate = 0.6)
	# shift times such that t=tstart --> t=1
	dat.no.trunc     <- subset(dat.full, tb >= tstart)
	dat.no.trunc$tb  <- dat.no.trunc$tb - tstart + 1
	# Truncates:
	dat.chopped <- subset(dat.no.trunc, tb <= ttrunc)
	
	# Set parameters for every models:
	PRM <- create.model.prm(dat      = dat.chopped,
							dat.full = dat.no.trunc,
							horizon  = horizon,
							horiz.fcast = horiz.fcast,  
							GI.mean = as.numeric(GI.mean), 
							GI.stdv = as.numeric(GI.stdv),
							GI.dist =  "gamma",
							cori.window = cori.window,
							GI_span = GI_span,
							pop_size = pop_size,
							mcmc_iter = mcmc_iter,
							mcmc_nchains = mcmc_nchains,
							mcmc_diagnostic = mcmc_diagnostic,
							SEmInR.prm.fxd = SEmInR.prm.fxd,
							SEmInR.prm.to.fit = SEmInR.prm.to.fit)
	
	# Forecast:
	
	# If just one forecasting model (typically used for test/debug):
	if(!is.null(single.model.fcast)) {
		PRM <- PRM[[single.model.fcast]]
		fcast <- fcast_incidence(prms = PRM, do.plot = do.plot)
	}
	if(is.null(single.model.fcast)) {
		fcast <- try(lapply(X = PRM,
							FUN = fcast_incidence,
							do.plot = do.plot),
					 silent = TRUE)
		# # # DEBUG: removed 'try'
		# fcast <- lapply(X = PRM,
		# 				FUN = fcast_incidence,
		# 				do.plot = do.plot)
		# - - - - 
	}
	return(fcast)
}


get.synthetic.data.db <- function(db.path,
								  country = NULL,
								  disease = NULL,
								  synthetic = NULL,
								  source.keys = NULL,
								  eventtype = NULL,
								  eventtype2 = NULL,
								  social.struct = NULL,
								  save.to.file = TRUE){
	# Get the incidence curves from data base
	# (there are potentially several Monte-Carlo realizations)
	
	# Load synthetic data:
	dat0 <- read.database(db.path,
						  country,
						  disease,
						  synthetic,
						  source.keys,
						  eventtype,
						  eventtype2,
						  social.struct)
	# Reformat (keeps only what's used later on)
	inc.tb <- convert.for.backtest(dat0)
	return(inc.tb)
}


get.synthetic.epi.param <- function(source){
	# Extract parameter information 
	# located in column named 'source'
	prm.name <- get.prm.names.from.source(source.string = source)
	x <- vector()
	for(i in 1:length(prm.name)){
		x[i] <- get.prm.value.from.source(source, prm.name = prm.name[i])[1]
	}
	names(x) <- prm.name
	return(x)
}


get.GI <- function(trueparam){
	### Retrieve GI from model parameters
	###
	seir <- sum(grepl('DOL.days',names(trueparam)))
	if(seir>0) {
		DOL.true <- trueparam['DOL.days']
		DOI.true <- trueparam['DOI.days']
		nI.true <- trueparam['nI']
		nE.true <- trueparam['nI']
		GI.mean <- (DOL.true+DOI.true/2/nI.true*(nI.true+1))
		GI.var <- GI.mean/mean(nI.true,nE.true) # <-- not sure about that!!!
	}
	if(seir==0) {
		GI.mean <- trueparam['GImean']
		GI.var <- trueparam['GIvar']
	}
	res <- as.numeric(c(GI.mean,GI.var))
	return(c(GI.mean=res[1], GI.var=res[2]))
}

get.trunc.time <- function(file,trueparam){
	
	# Parameters for the backtesting:
	prm <- read.csv(file,header = FALSE)
	
	# Truncation date (synthetic beyond
	# this time is supposed unknown)
	trunc.type <- as.character(prm[prm[,1]=="trunc.type",2])
	trunc.date <- as.numeric(as.character(prm[prm[,1]=="trunc.date",2]))
	trunc.gen  <- as.numeric(as.character(prm[prm[,1]=="trunc.gen",2]))
	if(trunc.type=="date") trunc.gen <- NULL
	if(trunc.type=="generation") trunc.date <- NULL
	if(trunc.type!="generation" & trunc.type!="date") trunc.gen <- trunc.date <- NULL
	
	# # Generation interval
	# seir <- sum(grepl('DOL.days',names(trueparam)))
	# if(seir>0) {
	# 	DOL.true <- trueparam['DOL.days']
	# 	DOI.true <- trueparam['DOI.days']
	# 	nI.true <- trueparam['nI']
	# 	GI.mean <- (DOL.true+DOI.true/2/nI.true*(nI.true+1))
	# }
	# if(seir==0)	 GI.mean <- trueparam['GImean']
	GI.mean <- get.GI(trueparam)['GI.mean']
	# truncated time:
	if(is.null(trunc.date)) ttr <- as.numeric(GI.mean*trunc.gen)
	if(is.null(trunc.gen))  ttr <- as.numeric(trunc.date)
	return(ttr)
}




read_prm <- function(file,x){
	f <- read.csv(file,header = F)
	as.numeric(as.character(f[trimws(as.character(f[,1]))==x,2]))
}

read_all_prm <-function(){
	# Backtesting Parameters 
	fpb            <<- 'prm_backtest.csv'
	horizon        <<- read_prm(fpb,'horizon') 
	horiz.fcast    <<- read_prm(fpb,'horiz.fcast') 
	GI.bias        <<- read_prm(fpb,'GI.bias') 
	n.MC.max       <<- read_prm(fpb,'n.MC.max')
	CI             <<- read_prm(fpb,'CI') 
	multicores     <<- read_prm(fpb,'parallel')
	rel.err        <<- read_prm(fpb,'relError')
	cori.window    <<- read_prm(fpb,'cori.window') # <-- TO DO: should be in another file!
	
	# RESuDe parameters:
	fresude         <<- 'prm-resude.csv'
	GI_span         <<- read_prm(file = fresude, x='GI_span')
	pop_size        <<- read_prm(fresude,'pop_size')
	mcmc_iter       <<- read_prm(fresude,'mcmc_iter')
	mcmc_nchains    <<- read_prm(fresude,'mcmc_nchains')
	mcmc_diagnostic <<- read_prm(fresude,'mcmc_diagnostic')
	
	# SEmInR parameters:
	fseminr      <<- 'prm-seminr.csv'
	horizon      <<- read_prm(fseminr, 'horizon')
	nE           <<- read_prm(fseminr, 'nE')
	nI           <<- read_prm(fseminr, 'nI')
	init_I1      <<- read_prm(fseminr, 'init_I1')
	n.time.steps <<- read_prm(fseminr, 'n.time.steps')
	per.capita   <<- read_prm(fseminr, 'per.capita')
	
	latent_mean       <<- read_prm(fseminr, 'latent_mean')
	infectious_mean   <<- read_prm(fseminr, 'infectious_mean')
	popSize           <<- read_prm(fseminr, 'popSize')
	R0                <<- read_prm(fseminr, 'R0')
	
	SEmInR.prm.to.fit <<- c(
		latent_mean     = latent_mean,
		infectious_mean = infectious_mean,
		popSize         = popSize,
		R0              = R0)
	
	SEmInR.prm.fxd <<-  c(
		horizon = horizon,
		nE = nE,
		nI = nI,
		init_I1 = init_I1,
		n.time.steps = n.time.steps,
		per.capita = per.capita
	)	
}
