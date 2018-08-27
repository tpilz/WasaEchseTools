#' Wrapper function for calibration of the WASA-SED model
#'
#' A wrapper function, executing a simulation of the WASA-SED model for a given
#' parameter realisation.
#'
#' @param pars A named vector of type numeric with the values of selected parameters
#' (directed to \code{\link[WasaEchseTools]{wasa_modify_pars}}).
#'
#' @param wasa_app Character string giving the system command of the WASA-SED application
#' (directed to \code{\link[WasaEchseTools]{wasa_run}}).
#'
#' @param sp_input_dir Character string of the directory containing the spatial WASA-SED
#' input. E.g. the output of \code{\link[lumpR]{db_wasa_input}} (argument 'dest_dir')
#' (directed to \code{\link[WasaEchseTools]{wasa_prep_runs}}).
#'
#' @param meteo_dir Character string of the directory containing the time series input
#' files for the WASA-SED model. Requires the files temperature.dat, radiation.dat,
#' humidity.dat, and rain_daily.dat or rain_hourly.dat (directed to \code{\link[WasaEchseTools]{wasa_prep_runs}}).
#'
#' @param dir_run Character specifying the directory for the model run with the current
#' parameter realisation (directed to \code{\link[WasaEchseTools]{wasa_run}}).
#' Default: A temporary directory created with \code{\link{tempfile}}.
#'
#' @param sim_start Object of class 'date' giving the start date of the simulation
#' (will be written into WASA-SED input file 'do.dat'; directed to \code{\link[WasaEchseTools]{wasa_prep_runs}}).
#'
#' @param sim_end Object of class 'date' giving the end date of the simulation
#' (will be written into WASA-SED input file 'do.dat' directed to \code{\link[WasaEchseTools]{wasa_prep_runs}}).
#'
#' @param resol Temporal resolution of the model simulations. Supported are: 'daily'
#' (default) and 'hourly'.
#'
#' @param warmup_start An object of class 'date' giving the start date of the warm-up period.
#' If \code{NULL} (default), the value in do.dat (lines 4 and 6) is used.
#' Directed to \code{\link[WasaEchseTools]{wasa_run}}.
#'
#' @param radex_file Character string of the file (including path) of preprared
#' extraterrestrial radiation input (named 'extraterrestrial_radiation.dat')
#' (directed to \code{\link[WasaEchseTools]{wasa_prep_runs}})).
#' 
#' @param dat_streamflow OPTIONAL: Object of class 'xts' containing a time series of
#' streamflow at the catchment outlet in (m3/s) in the resolution of the model run.
#' NA values are allowed and will be discarded for goodness of fit calculations.
#' Only needed if \code{return_val = 'nse'}.
#'
#' @param dat_pr OPTIONAL: Object of class 'xts' containing a time series of catchment-wide
#' average precipitation in (m3) in the resolution of the model run. If not given (default),
#' it will be read (and calculated) from the WASA-SED input files if needed. However,
#' it is more efficient in terms of function execution time to specify it as input if needed.
#' Only needed if \code{return_val = 'hydInd'}.
#'
#' @param flood_thresh OPTIONAL: Numeric value giving the threshold in (m3/s) for the definition of
#' a flood event (directed to \code{\link[WasaEchseTools]{hydInd}}). Only needed if
#' \code{return_val = 'hydInd'}.
#'
#' @param thresh_zero OPTIONAL: Values of discharge in (m3/s) below this value will be treated as
#' zero flows (directed to \code{\link[WasaEchseTools]{hydInd}}). Only needed if
#' \code{return_val = 'hydInd'}.
#'
#' @param warmup_len Integer giving the length of the warm-up period in months. Default: 3.
#' Directed to \code{\link[WasaEchseTools]{wasa_run}}.
#'
#' @param max_pre_runs Integer specifying the maximum number of warm-up iterations to be
#' applied. If the relative storage change is still larger than \code{storage_tolerance}
#' after \code{max_pre_runs} iterations, the warm-up will be aborted and the model be run
#' anyway. A warning will be issued. Default: 20. Directed to \code{\link[WasaEchseTools]{wasa_run}}.
#'
#' @param storage_tolerance Numeric value giving the relative change of the model's water
#' storages between two connsecutive warm-up runs below which the warm-up will be
#' concluded and the actual model simulation be started. Default: 0.01.
#' Directed to \code{\link[WasaEchseTools]{wasa_run}}.
#'
#' @param return_val Character string specifying your choice of what this function
#' shall return. Default: 'river_flow'. See description of return value below.
#'
#' @param keep_rundir Value of type \code{logical}. Shall directory \code{dir_run}
#' be retained (\code{TRUE}) or deleted (\code{FALSE}) after function execution?
#' Default: \code{FALSE}.
#'
#' @param keep_log Value of type \code{logical}. Shall a log file of the model run be written (to \code{dir_run})?
#' Default: \code{FALSE}. Will be ignored if \code{keep_rundir = FALSE}.
#' Directed to \code{\link[WasaEchseTools]{wasa_run}}.
#'
#' @details Function is a wrapper function, internally executing functions \code{\link[WasaEchseTools]{wasa_prep_runs}},
#' \code{\link[WasaEchseTools]{wasa_modify_pars}}, and \code{\link[WasaEchseTools]{wasa_run}}
#' and processing and returning the simulation output as specified. The function can
#' be employed by model calibration functions such as \code{\link[HydroBayes]{dream}}
#' or \code{\link[ppso]{optim_dds}}, or to execute single model runs within a single
#' function call.
#'
#' @return Function returns a vector of numeric values. Can be controlled by argument
#' \code{return_val}. Currently implemented are the options:
#'
#' river_flow: An object of class 'xts' containing the simulated river flow leaving the
#' catchment outlet in m3/s for the specified simulation period and resolution.
#'
#' hydInd: Named vector of hydrological indices calculated with function \code{\link[WasaEchseTools]{hydInd}}.
#' See function's doc for more information. This option requires the optional input arguments
#' \code{flood_thresh}, \code{thresh_zero}, and (optional) \code{dat_pr}.
#' 
#' nse: A single numeric value giving the Nash-Sutcliffe efficiency of streamflow
#' simulation at the catchment outlet (a common goodness of fit measure of hydrological
#' model runs).
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
wasa_calibwrap <- function(
  pars = NULL,
  wasa_app = NULL,
  sp_input_dir = NULL,
  meteo_dir = NULL,
  dir_run = paste0(tempfile(pattern = "wasa_calib_"), sep="/"),
  sim_start = NULL,
  sim_end = NULL,
  resol = "daily",
  warmup_start = NULL,
  radex_file = NULL,
  dat_streamflow = NULL,
  dat_pr = NULL,
  flood_thresh = NULL,
  thresh_zero = NULL,
  warmup_len = 3,
  max_pre_runs = 20,
  storage_tolerance = 0.01,
  return_val = "river_flow",
  keep_rundir = FALSE,
  keep_log = FALSE
) {
  # CHECKS #
  if(resol == "daily") {
    timestep  <- 24
    prec_file <- "rain_daily.dat"
  } else if(resol == "hourly") {
    timestep <- 1
    prec_file <- "rain_hourly.dat"
  }
  if(return_val == "hydInd") {
    if(!is.numeric(flood_thresh)) stop("Argument return_val = hydInd requires argument flood_thresh to be given!")
    if(!is.numeric(thresh_zero)) stop("Argument return_val = hydInd requires argument thresh_zero to be given!")
  }
  if(return_val == "nse") {
    if(!is.xts(dat_streamflow)) stop("Argument return_val = nse requires an object of class 'xts' for argument dat_streamflow!")
  }

  # prepare run directory
  wasa_prep_runs(sp_input_dir = sp_input_dir, meteo_dir = meteo_dir,
                 radex_file = radex_file, wasa_sim_dir = dir_run, proj_name = "",
                 sim_start = sim_start, sim_end = sim_end, timestep = timestep,
                 outfiles ="river_flow")

  # modify wasa input
  wasa_modify_pars(pars, paste(dir_run, "input", sep="/"))

  # run wasa (including warmup)
  wasa_run(dir_run, wasa_app, warmup_start, warmup_len, max_pre_runs, storage_tolerance, keep_log = keep_log)

  # get simulations
  file_wasa <- paste(dir_run, "output/River_Flow.out", sep="/")
  dat_wasa <- read.table(file_wasa, header=T, skip=1, check.names = F) %>%
    mutate(date = as.POSIXct(paste(year, day, sep="-"), "%Y-%j", tz ="UTC"), group = "wasa") %>%
    rename(value = "1") %>%
    dplyr::select(date, group, value)

  if(return_val == "hydInd") {
    dat_sim_xts <- xts(dat_wasa$value, dat_wasa$date)
    # precipitation (model forcing); get catchment-wide value (area-weighted precipitation mean)
    if(is.null(dat_pr)) {
      dat_sub <- readLines(paste(dir_run, "input/Hillslope/hymo.dat", sep="/"))
      dat_sub <- dat_sub[-c(1,2)]
      dat_sub_area <- sapply(dat_sub, function(x) as.numeric(unlist(strsplit(x, "\t"))[c(1,2)]), USE.NAMES = F)
      dat_sub_area <- dat_sub_area[2, order(dat_sub_area[1,])]
      sub_pars <- data.frame(object = paste0("sub_", 1:length(dat_sub_area)), area = dat_sub_area)
      dat_prec <- left_join(read.table(paste(meteo_dir, prec_file, sep="/"), header=T, skip=2, check.names = F)[,-2] %>%
                              mutate(date = as.POSIXct(sprintf(.[[1]], fmt="%08d"), "%d%m%Y", tz="UTC")) %>%
                              dplyr::select(-1) %>%
                              melt(id.vars="date", variable.name = "object") %>%
                              mutate(object = paste0("sub_", object), variable = "precip"),
                            sub_pars %>%
                              mutate_if(is.factor, as.character),
                            by = "object") %>%
        # in case of hourly resolution, calculate daily sums (WASA output River_Flow.out is only daily, even in case of hourly simulation resolution)
        group_by(date, variable, object, area) %>%
        summarise(value = sum(value)) %>%
        group_by(date, variable) %>%
        # sums over the catchment in m3
        summarise(value = sum(value * area*1e3)) %>%
        filter(date >= as.POSIXct(paste(sim_start, "00:00:00")) & date <= as.POSIXct(paste(sim_end, "23:59:59")))
      dat_pr <- xts(dat_prec$value, dat_prec$date)
    }
    # calculate diagnostic values
    out_vals <- suppressWarnings(hydInd(dat_sim_xts, dat_pr, na.rm = T, thresh.zero = thresh_zero, flood.thresh = flood_thresh))
    out_vals[which(is.na(out_vals) | is.nan(out_vals))] <- 0
    
  } else if(return_val == "river_flow") {
    out_vals <- xts(dat_wasa$value, dat_wasa$date)
    
  } else if(return_val == "nse") {
    dat_nse <- left_join(select(dat_wasa, date, value),
                         data.frame(obs = dat_streamflow,
                                    date = index(dat_streamflow)),
                         by = "date") %>%
      mutate(diffsq = (value - obs)^2,
             obsmean = mean(obs), diffsqobs = (obs - obsmean)^2) %>%
      summarise(nse = 1 - sum(diffsq) / sum(diffsqobs))
    out_vals <- as.numeric(dat_nse)
    
  }

  # clean up
  if(!keep_rundir) unlink(dir_run, recursive = T)

  # output
  return(out_vals)
} # EOF
