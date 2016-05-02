# backtest_fcast
Backtesting incidence forecasting of various models


# WARNINGS

This repo assumes the following repos are cloned in sister directories as this one:
 
* [Datsid](https://github.com/Outbreak-analysis/Datsid): Database of infectious series time series 
* [RESuDe](https://github.com/davidchampredon/RESuDe_forecast): Model implementing a renewal equation with susceptible depletion

This repo is still a work in progress... contact David Champredon for any questions.


# Execution

If the `Datsid` database is empty (typically when the repo has just been cloned): use the script 'runall'. It will first populate the database named `bcktst.db` with synthetic data, then perform the backtesting of forecast. 
The command line is: `runall n`, with `n` the number of Monte-Carlo iteration for the synthetic data.

It the `Datsid/bcktst.db` database is already populated with synthetic data, the command line is simply: `Rscript backtest_fcast.R`


 
# Code Structure

The main R script is `backtest_fcast.R` where the main function calling all the backtesting procedures is `fcast.wrap.mc`. 

To add (or remove) a forecasting model, edit the functions `create.model.prm` and `fcast_incidence`.