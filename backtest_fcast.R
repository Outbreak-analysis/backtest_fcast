###
###   FORECAST MODELS BACKTESTING
###

library(ggplot2);theme_set(theme_bw())
library(gridExtra)
library(plyr)
library(snowfall)
library(parallel)

t1 <- as.numeric(Sys.time())

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

# Models that generated synthetic data:
syn.models <- list("SEmInR", "RESuDe")

# Identify the source names of  synthetic data
db.path <- "../Datsid/bcktest.db"
bcktest <- get.list.sources(db.path = db.path)
idx <- lapply(syn.models, grepl, x = bcktest)
idx <- rowSums(matrix(unlist(idx),ncol=length(syn.models)))
# DEBUG bcktest <- bcktest[12:13]
bcktest <- bcktest[as.logical(idx)]
n.bcktest <- length(bcktest)

read_prm <- function(file,x){
	f <- read.csv(file,header = F)
	as.numeric(as.character(f[trimws(as.character(f[,1]))==x,2]))
}
# Backtesting Parameters 
fpb            <- 'prm_backtest.csv'
horizon        <- read_prm(fpb,'horizon') 
horiz.fcast    <- read_prm(fpb,'horiz.fcast') 
GI.bias        <- read_prm(fpb,'GI.bias') 
n.MC.max       <- read_prm(fpb,'n.MC.max')
CI             <- read_prm(fpb,'CI') 
multicores     <- read_prm(fpb,'parallel')
rel.err        <- read_prm(fpb,'relError')
cori.window    <- read_prm(fpb,'cori.window') # <-- TO DO: should be in another file!

# RESuDe parameters:
fresude         <- 'prm-resude.csv'
GI_span         <- read_prm(file = fresude, x='GI_span')
pop_size        <- read_prm(fresude,'pop_size')
mcmc_iter       <- read_prm(fresude,'mcmc_iter')
mcmc_nchains    <- read_prm(fresude,'mcmc_nchains')
mcmc_diagnostic <- read_prm(fresude,'mcmc_diagnostic')

# SEmInR parameters:
fseminr      <- 'prm-seminr.csv'
horizon      <- read_prm(fseminr, 'horizon')
nE           <- read_prm(fseminr, 'nE')
nI           <- read_prm(fseminr, 'nI')
init_I1      <- read_prm(fseminr, 'init_I1')
n.time.steps <- read_prm(fseminr, 'n.time.steps')
per.capita   <- read_prm(fseminr, 'per.capita')

latent_mean       <- read_prm(fseminr, 'latent_mean')
infectious_mean   <- read_prm(fseminr, 'infectious_mean')
popSize           <- read_prm(fseminr, 'popSize')
R0                <- read_prm(fseminr, 'R0')

SEmInR.prm.to.fit <- c(
	latent_mean     = latent_mean,
	infectious_mean = infectious_mean,
	popSize         = popSize,
	R0              = R0)

SEmInR.prm.fxd <-  c(
	horizon = horizon,
	nE = nE,
	nI = nI,
	init_I1 = init_I1,
	n.time.steps = n.time.steps,
	per.capita = per.capita
)


### 
### --- Run the backtesting ---
### 
sc.tmp <- list()
for(i in 1:n.bcktest){
	# Retrieve all synthetic epidemics from a model parameter set:
	dat.all <- get.synthetic.data.db(db.path     = db.path,
									 source.keys = bcktest[i],
									 eventtype   = 'incidence')
	mcvec <- unique(dat.all$mc)
	source <- dat.all$source[1]
	print(paste(i,"/",n.bcktest,":",source))
	# The 'true' parameters that generated these epidemics:
	trueparam  <- get.synthetic.epi.param(source = source)
	GI.mean    <- get.GI(trueparam)['GI.mean']
	GI.stdv    <- sqrt(get.GI(trueparam)['GI.var'])
	
	# Find time (after start date) 
	# to truncate full data (to make forecasts):
	ttrunc     <- ceiling(get.trunc.time(file = fpb,
									 trueparam = trueparam))
	# Parallel execution of the forecast
	# for a given scenario
	# across all MC realizations:
	n.cores <- detectCores()
	if(!multicores) n.cores <- 1
	sfInit(parallel = (n.cores>1), 
		   cpu = n.cores)
	sfLibrary(R0)
	sfLibrary(rstan)
	if(!use.DC.version.of.EpiEstim) sfLibrary(EpiEstim)
	
	idx.apply <- mcvec
	message(paste("Synthetic data contains",length(idx.apply),"MC iterations"))
	
	# Reduce backtesting to 
	# specified MC realizations:
	if(n.MC.max>0) {
		idx.apply <- idx.apply[1:n.MC.max]
		message(paste("but not more than",length(idx.apply),"are used."))
	}
	sfExportAll()
	res.parallel <- sfSapply(idx.apply, 
							 simplify = FALSE,
							 fcast.wrap.mc,
							 dat.all     = dat.all,
							 ttrunc      = ttrunc,
							 horizon     = horizon,
							 horiz.fcast = horiz.fcast,
							 GI.mean     = GI.mean,
							 GI.stdv     = GI.stdv,
							 cori.window = cori.window,
							 GI_span     = GI_span,
							 pop_size    = pop_size,
							 mcmc_iter   = mcmc_iter,
							 mcmc_nchains    = mcmc_nchains,
							 mcmc_diagnostic = mcmc_diagnostic,
							 SEmInR.prm.fxd    = SEmInR.prm.fxd,
							 SEmInR.prm.to.fit = SEmInR.prm.to.fit,
							 rel.err     = rel.err,
							 do.plot     = (n.cores==1)
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
# Merge and summarize scores for 
# all bactesting scenarios:
scsum <- merge.sum.scores(sc.tmp,CI)

# Plot:
plot.scores(scsum)

# Save all:
save.image('bcktst.RData')

# ==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
t2 <- as.numeric(Sys.time())
message(paste("Finished in",round((t2-t1)/60,2),"min"))
