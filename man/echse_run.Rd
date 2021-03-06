% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/echse_run.R
\name{echse_run}
\alias{echse_run}
\title{Executes simulations with the WASA engine of the ECHSE simulation environment}
\usage{
echse_run(dir_input = NULL, run_name = NULL, echse_app = NULL,
  choices = NULL, output_sel = NULL, output_dbg = NULL,
  sim_start = NULL, sim_end = NULL, resolution = 86400,
  warmup_start = NULL, warmup_len = 3, max_pre_runs = 20,
  storage_tolerance = 0.01, error2warn = FALSE, nthreads = 1)
}
\arguments{
\item{dir_input}{Character string giving the main simulation directory containing a
directory 'data' with the readily prepared ECHSE input (e.g. prepared with function
\code{\link[WasaEchseTools]{echse_prep_runs}}).}

\item{run_name}{Character string defining the name of a sub-directory within \code{dir_input}
to be created as storage of the actual simulations. It will further contain the directories
'pars' containing ECHSE's control files and parameter files with the model structure choices
(if \code{choices} is given), and 'out' containing the actual simulation output.}

\item{echse_app}{Character string giving the system command of the application.}

\item{choices}{A named data.frame with each element containing the flag for a specific
choice. See the latest version of ECHSE's WASA engine for required choice flags. If
\code{NULL} (the default), it is assumed the model is run with default or manually
defined selections. In that case, the 'data/parameter' input directory needs to
contain the actual sharedParamNum_WASA_* files instead of the *_tpl.dat files.}

\item{output_sel}{A named data.frame with elements 'object', 'variable', and 'digits'
defining the output to be generated. Needed for control file 'output_selection.txt'.
See the ECHSE core manual for further information. If \code{NULL} (default), the
variable 'out' of object 'node_su_out_1' will be given (i.e. the basin's outlet river
flow in m3/s).}

\item{output_dbg}{A character vector specifying the names of objects for which
debug output shall be written. Needed for control file 'output_debug.txt'. If set
to \code{NULL} (default), no debug output will be generated.}

\item{sim_start}{Character string giving the start date of the simulation period
in the format "\%Y-\%m-\%d \%H:\%M:\%S".}

\item{sim_end}{Character string giving the end date of the simulation period
in the format "\%Y-\%m-\%d \%H:\%M:\%S".}

\item{resolution}{Integer giving the simulation time step in seconds. Default:
86400 (i.e. daily resolution).}

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

\item{error2warn}{Value of type \code{logical}. Shall runtime errors of the model be
reported as a warning instead of stopping this function with an error? Default: \code{FALSE}.}

\item{nthreads}{Number of cores that shall be employed for the ECHSE run (argument
'number_of_threads' in configuration file). See ECHSE core manual for more information.}
}
\value{
Function returns nothing.
}
\description{
Function executes simulation with the WASA engine of the ECHSE simulation environment
including warm-up runs to start from stable model states as initial conditions.
}
\note{
To avoid warm-up runs, set \code{max_pre_runs} or \code{warmup_len} to zero.
}
\author{
Tobias Pilz \email{tpilz@uni-potsdam.de}
}
