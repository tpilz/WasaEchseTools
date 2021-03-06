% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/echse_prep_runs.R
\name{echse_prep_runs}
\alias{echse_prep_runs}
\title{Preparation of input for the WASA engine of the ECHSE simulation environment}
\usage{
echse_prep_runs(sp_input_dir = NULL, wgt_file = NULL,
  meteo_ext_datafiles = NULL, prep_tpl = FALSE,
  echse_sim_dir = getwd())
}
\arguments{
\item{sp_input_dir}{Character string of the directory containing the spatial input.
E.g. the output of \code{\link[lumpR]{db_echse_input}}.}

\item{wgt_file}{Character string specifying the name of the file containing the
weighting of input locations of meteorological data. Should be the output of function
\code{\link[geostat]{externalInputLocationsTable}} (argument 'file_result').}

\item{meteo_ext_datafiles}{Data.frame containing information of ECHSE's external
datafile input. Requires the variables 'variable', 'sums', 'past', and 'file'
(see the ECHSE manual for more information).}

\item{prep_tpl}{Logical. Shall template parameter files be created? If \code{TRUE},
files sharedParamNum_WASA_svc_tpl.dat, sharedParamNum_WASA_lu_tpl.dat, and
sharedParamNum_WASA_rch_tpl.dat will be created (and the old files be deleted
accordingly) where all choice_* parameter values will be replaced by placeholders
to be adjusted, for instance, within a multi-hypothesis study employing different
model structures (i.e. using different choice_* parameter realisations). Default: \code{FALSE}.}

\item{echse_sim_dir}{Character string of the location the output shall be written
to (directory will be created if it does not exist).}
}
\value{
Function returns nothing.
}
\description{
Makes use of spatial model input prepared with the \code{\link[lumpR]{lumpR}} package,
pre-processed meteorological data, and further input to this function to generate a
directory for the model ready to run simulations.
}
\author{
Tobias Pilz \email{tpilz@uni-potsdam.de}
}
