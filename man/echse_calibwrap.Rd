% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/echse_calibwrap.R
\name{echse_calibwrap}
\alias{echse_calibwrap}
\title{Wrapper function for calibration of the WASA engine of the ECHSE environment}
\usage{
echse_calibwrap(pars = NULL, dir_input = NULL,
  dir_run = paste0(tempfile(pattern = "echse_calib_"), sep = "/"),
  echse_app = NULL, choices = NULL, sim_start = NULL, sim_end = NULL,
  resolution = 86400, dat_pr = NULL, flood_thresh = NULL,
  thresh_zero = NULL, warmup_start = NULL, warmup_len = 3,
  max_pre_runs = 20, storage_tolerance = 0.01, return_val = "river_flow",
  keep_rundir = FALSE)
}
\arguments{
\item{pars}{A named vector of type numeric with the values of selected parameters.}

\item{dir_input}{Character string giving the main simulation directory containing a
directory 'data' with the readily prepared ECHSE input (e.g. prepared with function
\code{\link[WasaEchseTools]{echse_prep_runs}}).}

\item{dir_run}{Character specifying the directory for the model run with the current
parameter realisation. Default: A temporary directory created with \code{\link{tempfile}}.}

\item{echse_app}{Character string giving the system command of the application.}

\item{choices}{A named data.frame with each element containing the flag for a specific
choice. See the latest version of ECHSE's WASA engine for required choice flags. If
\code{NULL} (the default), it is assumed the model is run with default or manually
defined selections. In that case, the 'data/parameter' input directory needs to
contain the actual sharedParamNum_WASA_* files instead of the *_tpl.dat files.}

\item{sim_start}{Character string giving the start date of the simulation period
in the format "\%Y-\%m-\%d \%H:\%M:\%S".}

\item{sim_end}{Character string giving the end date of the simulation period
in the format "\%Y-\%m-\%d \%H:\%M:\%S".}

\item{resolution}{Integer giving the simulation time step in seconds. Default:
86400 (i.e. daily resolution).}

\item{dat_pr}{Object of class 'xts' containing a time series of catchment-wide
average precipitation in (m3) in the resolution of the model run. If not given (default),
it will be read (and calculated) from the ECHSE input files if needed. However,
it is more efficient in terms of function execution time to specify it as input if needed.
Only needed if \code{return_val = 'hydInd'}.}

\item{flood_thresh}{Numeric value giving the threshold in (m3/s) for the definition of
a flood event (directed to \code{\link[WasaEchseTools]{hydInd}}). Only needed if
\code{return_val = 'hydInd'}.}

\item{thresh_zero}{Values of discharge in (m3/s) below this value will be treated as
zero flows (directed to \code{\link[WasaEchseTools]{hydInd}}). Only needed if
\code{return_val = 'hydInd'}.}

\item{warmup_start}{Character string giving the start date of the warm-up period
in the format "\%Y-\%m-\%d \%H:\%M:\%S". If \code{NULL} (default), argument 'sim_start'
will be used.}

\item{warmup_len}{Integer giving the length of the warm-up period in months. Default: 3.}

\item{max_pre_runs}{Integer specifying the maximum number of warm-up iterations to be
applied. If the relative storage change is still larger than \code{storage_tolerance}
after \code{max_pre_runs} iterations, the warm-up will be aborted and the model be run
anyway. A warning will be issued. Default: 20.}

\item{storage_tolerance}{Numeric value giving the relative change of the model's water
storages between two connsecutive warm-up runs below which the warm-up will be
concluded and the actual model simulation be started. Default: 0.01.}

\item{return_val}{Character string specifying your choice of what this function
shall return. Default: 'river_flow'. See description of return value below.}

\item{keep_rundir}{Value of type \code{logical}. Shall directory \code{dir_run}
be retained (\code{TRUE}) or deleted (\code{FALSE}) after function execution?
Default: \code{FALSE}.}
}
\value{
Function returns a vector of numeric values. Can be controlled by argument
\code{return_val}. Currently implemented are the options:

river_flow: Simulated river flow leaving the catchment outlet in m3/s for the
specified simulation period and resolution

hydInd: Named vector of hydrological indices calculated with function \code{\link[WasaEchseTools]{hydInd}}.
See function's doc for more information. This option requires the optional input arguments
\code{flood_thresh}, \code{thresh_zero}, and (optional) \code{dat_pr}.
}
\description{
A wrapper function, executing a simulation of the WASA engine for a given
parameter realisation.
}
\details{
The function can be employed by model calibration functions such as
\code{\link[HydroBayes]{dream}} or \code{\link[ppso]{optim_dds}}, or to execute
single model runs within a single function call.
}
\author{
Tobias Pilz \email{tpilz@uni-potsdam.de}
}