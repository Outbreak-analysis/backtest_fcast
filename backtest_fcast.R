###
###   FORECAST MODELS BACKTESTING
###

t1 <- as.numeric(Sys.time())

library(ggplot2);theme_set(theme_bw())
library(gridExtra)
library(plyr)
library(snowfall)
library(parallel)

# - - - - - - - - - - - - - - - - - - - - - - -
# There is a problem in 'OverallInfectivity' 
# function when data set is short.
use.DC.version.of.EpiEstim <- TRUE  
if(!use.DC.version.of.EpiEstim) library(EpiEstim)
if(use.DC.version.of.EpiEstim) source("EstimationR.R")
# - - - - - - - - - - - - - - - - - - - - - - -
source("read-data.R")
source("forecast_fitsim.R")
source("forecast_utils.R")
source("scores.R")


# ==== DATA ====

# Models that generated synthetic data:
syn.models <- list("SEmInR", "RESuDe")

# Identify the source names of synthetic data:
db.path <- "../Datsid/bcktest.db"
bcktest <- get.list.sources(db.path = db.path)
idx     <- lapply(syn.models, grepl, x = bcktest)
idx     <- rowSums(matrix(unlist(idx),ncol=length(syn.models)))

# Subset the selected data sources:
bcktest <- bcktest[as.logical(idx)]
mydebug <- FALSE
if(mydebug)  bcktest <- bcktest[c(1,5,11,14)]
n.bcktest <- length(bcktest)


### ==== RUN BACKTESTING ====

sc.tmp <- list() # list for the scores

for(i in 1:n.bcktest){
	# Retrieve all synthetic epidemics from a model parameter set:
	dat.all <- get.synthetic.data.db(db.path     = db.path,
									 source.keys = bcktest[i],
									 eventtype   = 'incidence')
	mcvec  <- unique(dat.all$mc)
	source <- dat.all$source[1]
	message(paste(i,"/",n.bcktest,":",source))
	
	# The 'true' parameters that generated these epidemics:
	trueparam  <- get.synthetic.epi.param(source = source)
	GI.mean    <- get.GI(trueparam)['GI.mean']
	GI.stdv    <- sqrt(get.GI(trueparam)['GI.var'])
	
	# Read all parameters for backtest and forecasting models:
	read_technical_prm()
	set_true_param_for_not_fitted(trueparam)
	
	# Initialize starting values for models needing a minimization
	xinit<-vector()
	
	# Helping a bit to help convergence on the hundreds different data sets:
	xinit[['SEmInR_R0_init']]      <- runif(1,0.8,1.2) * trueparam[['R0']]
	xinit[['SEmInR_popsize_init']] <- runif(1,0.8,1.2) * trueparam[['pop.size']]
	xinit[['RESuDe_popsize_init']] <- runif(1,0.8,1.2) * trueparam[['pop.size']]
	init.prm.to.fit(xinit)
	
	# Find time (after start date) to truncate full data (to make forecasts):
	ttrunc <- ceiling(get.trunc.time(file = fpb,
									 trueparam = trueparam))
	
	# Parallel execution of the forecast for a given scenario
	# across all MC realizations:
	n.cores <- parallel::detectCores()
	if(!multicores) n.cores <- 1
	sfInit(parallel = (n.cores>1), 
		   cpu = n.cores)
	sfLibrary(R0)
	sfLibrary(rstan)
	sfLibrary(deSolve)
	if(!use.DC.version.of.EpiEstim) sfLibrary(EpiEstim)
	
	idx.apply <- mcvec
	message(paste("Synthetic data contains",length(idx.apply),"MC iterations"),appendLF = F)
	
	# Reduce backtesting to specified MC realizations (e.g., to save time):
	if(n.MC.max>0) {
		idx.apply <- idx.apply[1:n.MC.max]
		message(paste(" but not more than",length(idx.apply),"are used."))
	}
	sfExportAll()
	
	res.parallel <- sfSapply(x                 = idx.apply, 
							 simplify          = FALSE,
							 fun               = fcast.wrap.mc,
							 dat.all           = dat.all,
							 ttrunc            = ttrunc,
							 horizon           = horizon,
							 horiz.fcast       = horiz.fcast,
							 GI.mean           = GI.mean,
							 GI.stdv           = GI.stdv,
							 cori.window       = cori.window,
							 GI_span           = GI_span,
							 mcmc_iter         = mcmc_iter,
							 mcmc_nchains      = mcmc_nchains,
							 mcmc_diagnostic   = mcmc_diagnostic,
							 SEmInR.prm.fxd    = SEmInR.prm.fxd,
							 SEmInR.prm.to.fit = SEmInR.prm.to.fit,
							 rel.err           = rel.err,
							 do.plot           = (n.cores==1)
	)
	sfStop()

	### Calculate scores
	sc.tmp[[i]] <- calc.scores(res.parallel, horiz.fcast, rel.err)
	
	sc.tmp[[i]]$modelsyndata <- substr(x = source, start = 1, stop=6)
	sc.tmp[[i]]$source       <- source
	sc.tmp[[i]]$R0           <- trueparam['R0']
	sc.tmp[[i]]$GI.mean      <- GI.mean
	sc.tmp[[i]]$nMCdata      <- length(mcvec)
}

# Merge and summarize scores for all bactesting scenarios:
scsum <- merge.sum.scores(sc.tmp,CI)
plot.scores(scsum)

# Save all:
save.image('bcktst.RData')

# ==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
t2 <- as.numeric(Sys.time())
message(paste("Finished in",round((t2-t1)/60,2),"min"))
