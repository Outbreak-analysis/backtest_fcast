###
###   FORECAST _ONE_ SCENARIO WITH _ONE_ MODEL
###

library(ggplot2);theme_set(theme_bw())
library(gridExtra)
library(plyr)
library(snowfall)
library(parallel)

t1 <- as.numeric(Sys.time())
set.seed(1234)

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
syn.models <- list("SEmInR_4")#list("SEmInR", "RESuDe")
mc.choose <- 3
single.model.fcast <- 'GGM' #'RESuDe' # 'SEmInRdet'

# Identify the source names of  synthetic data
db.path <- "../Datsid/bcktest.db"
bcktest <- get.list.sources(db.path = db.path)
idx <- lapply(syn.models, grepl, x = bcktest)
idx <- rowSums(matrix(unlist(idx),ncol=length(syn.models)))
bcktest <- bcktest[as.logical(idx)]
n.bcktest <- length(bcktest)

# Read all parameters for backtest 
# and forecasting models
read_all_prm()
multicores <- 0

### 
### --- Run the backtesting ---
### 
sc.tmp <- list()
i <- 1
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
if(!multicores) n.cores <- 1

idx.apply <- mcvec
message(paste("Synthetic data contains",length(idx.apply),
			  "MC iterations but only one is used [",mc.choose,"]"))

if(mc.choose) idx.apply<- mcvec[mc.choose]

res.single <- fcast.wrap.mc(m = idx.apply, 
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
							do.plot     = TRUE,
							single.model.fcast = single.model.fcast
)

### Calculate scores
res <- list(res.single)
names(res) <- single.model.fcast

sc.tmp[[i]] <- calc.scores(res.parallel = list(res), 
						   horiz.fcast, rel.err)

sc.tmp[[i]]$modelsyndata <- substr(x = source, start = 1, stop=6)
sc.tmp[[i]]$source       <- source
sc.tmp[[i]]$R0           <- trueparam['R0']
sc.tmp[[i]]$GI.mean      <- GI.mean
sc.tmp[[i]]$nMCdata      <- length(mcvec)
# Merge and summarize scores for 
# all bactesting scenarios:
scsum <- merge.sum.scores(sc.tmp,CI)


# ==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
t2 <- as.numeric(Sys.time())
message(paste("Finished in",round((t2-t1)/60,2),"min"))
