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
#' @param return_val Character vector specifying your choice of what this function
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
#' @param error2warn Value of type \code{logical}. Shall runtime errors of the model be
#' reported as a warning instead of stopping this function with an error? If so, the
#' model run's log file will be saved and, in case of reasonable model output, this wrapper
#' function will proceed as usual. If no reasonable output could be found,
#' a value \code{NA} will be returned. Default: \code{FALSE}.
#' Also directed to \code{\link[WasaEchseTools]{wasa_run}}.
#'
#' @details Function is a wrapper function, internally executing functions \code{\link[WasaEchseTools]{wasa_prep_runs}},
#' \code{\link[WasaEchseTools]{wasa_modify_pars}}, and \code{\link[WasaEchseTools]{wasa_run}}
#' and processing and returning the simulation output as specified. The function can
#' be employed by model calibration functions such as \code{\link[HydroBayes]{dream}}
#' or \code{\link[ppso]{optim_dds}}, or to execute single model runs within a single
#' function call.
#'
#' @return Function returns a named list or a single element (if \code{length(return_val) = 1})
#' of numeric value(s). Can be controlled by argument \code{return_val}. Currently
#' implemented are the options:
#'
#' river_flow: An object of class 'xts' containing the simulated river flow leaving the
#' catchment outlet in m3/s for the specified simulation period and resolution.
#'
#' river_flow_mm: A single value of the catchment's simulated river outflow over
#' the specified simulation period in (mm).
#'
#' eta_mm: A single value of the catchment's simulated amount of actual evapotranspiration
#' over the specified simulation period in (mm).
#'
#' etp_mm: A single value of the catchment's simulated amount of potential evapotranspiration
#' over the specified simulation period in (mm).
#'
#' runoff_total_mm: A single value of the catchment's simulated amount of total runoff,
#' i.e. runoff contribution into river(s), over the specified simulation period in (mm).
#'
#' runoff_gw_mm: A single value of the catchment's simulated amount of groundwater runoff,
#' i.e. groundwater contribution into river(s), over the specified simulation period in (mm).
#'
#' runoff_sub_mm: A single value of the catchment's simulated amount of subsurface runoff,
#' i.e. near-surface soil water (excl. groundwater!) contribution into river(s), over the specified
#' simulation period in (mm). Note: due to computational reasons, 'runoff_gw_mm' will
#' be given in addition if this value is set.
#'
#' runoff_surf_mm: A single value of the catchment's simulated amount of surface runoff,
#' i.e. surface water contribution into river(s), over the specified simulation period in (mm).
#'
#' hydInd: Named vector of hydrological indices calculated with function \code{\link[WasaEchseTools]{hydInd}}.
#' See function's doc for more information. This option requires the optional input arguments
#' \code{flood_thresh}, \code{thresh_zero}, and (optional) \code{dat_pr}.
#'
#' nse: A single numeric value giving the Nash-Sutcliffe efficiency of streamflow
#' simulation at the catchment outlet (a common goodness of fit measure of hydrological
#' model runs).
#'
#' kge: A single numeric value giving the Kling-Gupta efficiency of streamflow
#' simulation at the catchment outlet (see paper \url{https://doi.org/10.1016/j.jhydrol.2009.08.003}).
#'
#' @note To avoid warm-up runs, set \code{max_pre_runs} or \code{warmup_len} to zero.
#'
#' Model's water balance: runoff_total_mm = runoff_surf_mm + runoff_sub_mm + runoff_gw_mm.
#' Moreover, river_flow_mm should be slight less than (due to transmission losses and storage
#' changes; but all in all almost equal to) runoff_total_mm as long as there is no
#' artificial influence (e.g. from reservoirs).
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
  keep_log = FALSE,
  error2warn = FALSE
) {
  # CHECKS #
  if(resol == "daily") {
    timestep  <- 24
    prec_file <- "rain_daily.dat"
  } else if(resol == "hourly") {
    timestep <- 1
    prec_file <- "rain_hourly.dat"
  }
  if("hydInd" %in% return_val) {
    if(!is.numeric(flood_thresh)) stop("Argument return_val = hydInd requires argument flood_thresh to be given!")
    if(!is.numeric(thresh_zero)) stop("Argument return_val = hydInd requires argument thresh_zero to be given!")
  }
  if(any(c("nse", "kge") %in% return_val)) {
    if(!is.xts(dat_streamflow)) stop("Argument return_val = nse requires an object of class 'xts' for argument dat_streamflow!")
  }

  # determine output files to be produced by WASA
  outfiles <- NULL; return_val_gr <- NULL
  if(any(c("nse", "kge", "hydInd", "river_flow", "river_flow_mm") %in% return_val))
    outfiles <- c(outfiles, "river_flow")

  if("runoff_surf_mm" %in% return_val) outfiles <- c(outfiles, "daily_total_overlandflow")
  if("runoff_sub_mm" %in% return_val) {
    outfiles <- c(outfiles, "daily_subsurface_runoff")
    # needed here as in WASA subsurface runoff = lateral near-surface soil water runoff + groundwater runoff!
    if(!("runoff_gw_mm" %in% return_val)) return_val <- c(return_val, "runoff_gw_mm")
  }
  if("runoff_gw_mm" %in% return_val) outfiles <- c(outfiles, "gw_discharge")

  if(any(c("eta_mm", "etp_mm", "runoff_total_mm", "runoff_sub_mm", "runoff_surf_mm") %in% return_val) && timestep < 24)
    stop("Currently a 'return_val' of 'eta_mm', 'etp_mm', or 'runoff_total_mm' can only be returned when running with daily resolution!")
  if("eta_mm" %in% return_val) outfiles <- c(outfiles, "daily_actetranspiration")
  return_val_gr <- c(return_val_gr, "daily_actetranspiration" = "eta_mm")
  if("etp_mm" %in% return_val) outfiles <- c(outfiles, "daily_potetranspiration")
  return_val_gr <- c(return_val_gr, "daily_potetranspiration" = "etp_mm")
  if("runoff_total_mm" %in% return_val) outfiles <- c(outfiles, "daily_water_subbasin")
  return_val_gr <- c(return_val_gr, "daily_water_subbasin" = "runoff_total_mm")

  # prepare run directory
  wasa_prep_runs(sp_input_dir = sp_input_dir, meteo_dir = meteo_dir,
                 radex_file = radex_file, wasa_sim_dir = dir_run, proj_name = "",
                 sim_start = sim_start, sim_end = sim_end, timestep = timestep,
                 outfiles = outfiles)

  # modify wasa input
  wasa_modify_pars(pars, paste(dir_run, "input", sep="/"))

  # run wasa (including warmup)
  wasa_run(dir_run, wasa_app, warmup_start, warmup_len, max_pre_runs, storage_tolerance,
           keep_log = keep_log, error2warn = error2warn)

  # get simulation data
  dat_wasa <- NULL
  if("river_flow" %in% outfiles) {
    file_wasa <- paste(dir_run, "output/River_Flow.out", sep="/")
    tryCatch(dat_file <- read.table(file_wasa, header=T, skip=1, check.names = F),
             error = function(e) stop(paste("Error reading", file_wasa, ":", e)))
    dat_wasa <- dat_file %>%
      mutate(date = as.POSIXct(paste(year, day, sep="-"), "%Y-%j", tz ="UTC"), group = "river_flow", subbas = -999) %>%
      rename(value = "1") %>%
      dplyr::select(date, group, value, subbas) %>%
      bind_rows(dat_wasa)
  }
  if("gw_discharge" %in% outfiles) {
    file_wasa <- paste(dir_run, "output/gw_discharge.out", sep="/")
    tryCatch(dat_file <- read.table(file_wasa, header=T, skip=1, check.names = F),
             error = function(e) stop(paste("Error reading", file_wasa, ":", e)))
    dat_wasa <- dat_file %>%
      select(Day, Year, num_range(prefix = "", range = 1:999)) %>%
      mutate(date = as.POSIXct(paste(Year, Day, sep="-"), "%Y-%j", tz ="UTC"), group = "gw_discharge") %>%
      select(-Day, -Year) %>%
      gather(key = subbas, value = value, -group, -date) %>%
      mutate(subbas = as.integer(subbas)) %>%
      bind_rows(dat_wasa)
  }
  if("daily_subsurface_runoff" %in% outfiles) {
    file_wasa <- paste(dir_run, "output/daily_subsurface_runoff.out", sep="/")
    tryCatch(dat_file <- read.table(file_wasa, header=T, skip=1, check.names = F),
             error = function(e) stop(paste("Error reading", file_wasa, ":", e)))
    dat_wasa <- dat_file %>%
      select(Day, Year, num_range(prefix = "", range = 1:999)) %>%
      mutate(date = as.POSIXct(paste(Year, Day, sep="-"), "%Y-%j", tz ="UTC"), group = "runoff_sub") %>%
      select(-Day, -Year) %>%
      gather(key = subbas, value = value, -group, -date) %>%
      mutate(subbas = as.integer(subbas)) %>%
      bind_rows(dat_wasa)
  }
  if("daily_total_overlandflow" %in% outfiles) {
    file_wasa <- paste(dir_run, "output/daily_total_overlandflow.out", sep="/")
    tryCatch(dat_file <- read.table(file_wasa, header=T, skip=1, check.names = F),
             error = function(e) stop(paste("Error reading", file_wasa, ":", e)))
    dat_wasa <- dat_file %>%
      select(Day, Year, num_range(prefix = "", range = 1:999)) %>%
      mutate(date = as.POSIXct(paste(Year, Day, sep="-"), "%Y-%j", tz ="UTC"), group = "runoff_surf") %>%
      select(-Day, -Year) %>%
      gather(key = subbas, value = value, -group, -date) %>%
      mutate(subbas = as.integer(subbas)) %>%
      bind_rows(dat_wasa)
  }
  if(any(c("daily_actetranspiration", "daily_potetranspiration", "daily_water_subbasin") %in% outfiles)) {
    outfiles_t <- outfiles[match(c("daily_actetranspiration", "daily_potetranspiration", "daily_water_subbasin"), outfiles, nomatch = 0)]
    dat_wasa <- map_dfr(outfiles_t, function(f) {
      file_wasa <- paste0(dir_run, "/output/", f, ".out")
      tryCatch(dat_file <- read.table(file_wasa, header=T, skip=1, check.names = F),
               error = function(e) stop(paste("Error reading", file_wasa, ":", e)))
      dat_file %>%
        mutate(date = as.POSIXct(paste(Year, Day, sep="-"), "%Y-%j", tz ="UTC"), group = return_val_gr[grep(f, names(return_val_gr))]) %>%
        select(-Day, -Year) %>%
        gather(key = subbas, value = value, -group, -date) %>%
        mutate(subbas = as.integer(subbas))
    }) %>%
      bind_rows(dat_wasa)
  }

  # ignore errors if desired
  if(nrow(dat_wasa) == 0) {
    if(!error2warn) {
      stop(paste0("There could be no output extracted from a model run, dir_run = ", dir_run))
    } else {
      warning(paste0("No reasonable model output could be extracted from model run in ", dir_run, ". NA will be returned!"))
      return(NA)
    }
  }

  # prepare output according to specifications
  out <- NULL
  if("hydInd" %in% return_val) {
    dat_tmp <- filter(dat_wasa, group == "river_flow")
    dat_sim_xts <- xts(dat_tmp$value, dat_tmp$date)
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

    out[["hydInd"]] <- out_vals

  }
  if("river_flow" %in% return_val) {
    dat_tmp <- filter(dat_wasa, group == "river_flow")
    out_vals <- xts(dat_tmp$value, dat_tmp$date)

    out[["river_flow"]] <- out_vals

  }
  if("river_flow_mm" %in% return_val) {
    dat_tmp <- filter(dat_wasa, group == "river_flow")
    # get catchment area (m2)
    dat_sub <- readLines(paste(dir_run, "input/Hillslope/hymo.dat", sep="/"))
    dat_sub <- dat_sub[-c(1,2)]
    dat_sub_area <- sapply(dat_sub, function(x) as.numeric(unlist(strsplit(x, "\t"))[c(1,2)]), USE.NAMES = F)
    dat_sub_area <- dat_sub_area[2, order(dat_sub_area[1,])]
    dat_sub_area <- sum(dat_sub_area) * 1e6
    # sum of catchment river outflow (mm)
    outflow <- dat_tmp$value * timestep * 60*60 # m3/day
    outflow <- sum(outflow) # m3
    outflow <- outflow *1000 # L
    outflow <- outflow / dat_sub_area # L/m2 = mm

    out[["river_flow_mm"]] <- outflow

  }
  if("nse" %in% return_val) {
    dat_tmp <- filter(dat_wasa, group == "river_flow")
    dat_nse <- left_join(select(dat_tmp, date, value),
                         data.frame(obs = dat_streamflow,
                                    date = index(dat_streamflow)),
                         by = "date") %>%
      drop_na(obs) %>%
      mutate(diffsq = (value - obs)^2,
             obsmean = mean(obs), diffsqobs = (obs - obsmean)^2) %>%
      summarise(nse = 1 - sum(diffsq) / sum(diffsqobs))
    out_vals <- as.numeric(dat_nse)

    out[["nse"]] <- out_vals

  }
  if("kge" %in% return_val) {
    dat_tmp <- filter(dat_wasa, group == "river_flow")
    dat_kge <- left_join(select(dat_tmp, date, value),
                         data.frame(obs = dat_streamflow,
                                    date = index(dat_streamflow)),
                         by = "date") %>%
      drop_na(obs) %>%
      summarise(a = cor(obs, value) - 1, b = sd(value)/sd(obs) - 1, c = mean(value)/mean(obs) - 1,
                kge = 1 - sqrt(a^2 + b^2 + c^2))
    out_vals <- as.numeric(dat_kge$kge)

    out[["kge"]] <- out_vals

  }
  if("runoff_gw_mm" %in% return_val) {
    dat_tmp <- filter(dat_wasa, group == "gw_discharge")
    # get catchment area (m2)
    dat_sub <- readLines(paste(dir_run, "input/Hillslope/hymo.dat", sep="/"))
    dat_sub <- dat_sub[-c(1,2)]
    dat_sub_area <- t(sapply(dat_sub, function(x) as.numeric(unlist(strsplit(x, "\t"))[c(1,2)]), USE.NAMES = F))
    colnames(dat_sub_area) <- c("subbas", "area")
    # merge data and calculate area-weighted mean and time period sum
    dat_gw_mm <- left_join(dat_tmp,
                          as.data.frame(dat_sub_area) %>%
                            mutate(subbas = as.integer(subbas), area_sum = sum(area)),
                          by = "subbas") %>%
      group_by(date, area_sum) %>%
      summarise(value_sum = sum(value)) %>% # sum over catchment (m3/timestep)
      mutate(value_mm = value_sum*1000 / (area_sum*1e6)) %>% # daily sum over catchment in (mm)
      ungroup() %>%
      summarise(value_sum_mm = sum(value_mm)) # sum over simulation period (mm)

    out[["runoff_gw_mm"]] <- dat_gw_mm$value_sum_mm

  }
  if("runoff_sub_mm" %in% return_val) {
    dat_tmp <- filter(dat_wasa, group == "runoff_sub")
    # get catchment area (m2)
    dat_sub <- readLines(paste(dir_run, "input/Hillslope/hymo.dat", sep="/"))
    dat_sub <- dat_sub[-c(1,2)]
    dat_sub_area <- t(sapply(dat_sub, function(x) as.numeric(unlist(strsplit(x, "\t"))[c(1,2)]), USE.NAMES = F))
    colnames(dat_sub_area) <- c("subbas", "area")
    # merge data and calculate area-weighted mean and time period sum
    dat_sub_mm <- left_join(dat_tmp,
                           as.data.frame(dat_sub_area) %>%
                             mutate(subbas = as.integer(subbas), area_sum = sum(area)),
                           by = "subbas") %>%
      group_by(date, area_sum) %>%
      summarise(value_sum = sum(value)) %>% # sum over catchment (m3/day)
      mutate(value_mm = value_sum*1000 / (area_sum*1e6)) %>% # daily sum over catchment in (mm)
      ungroup() %>%
      summarise(value_sum_mm = sum(value_mm)) # sum over simulation period (mm)

    out[["runoff_sub_mm"]] <- dat_sub_mm$value_sum_mm - out[["runoff_gw_mm"]] # NOTE: in WASA subsurface runoff = lateral near-surface soil water runoff + groundwater runoff!

  }
  if("runoff_surf_mm" %in% return_val) {
    dat_tmp <- filter(dat_wasa, group == "runoff_surf")
    # get catchment area (m2)
    dat_sub <- readLines(paste(dir_run, "input/Hillslope/hymo.dat", sep="/"))
    dat_sub <- dat_sub[-c(1,2)]
    dat_sub_area <- t(sapply(dat_sub, function(x) as.numeric(unlist(strsplit(x, "\t"))[c(1,2)]), USE.NAMES = F))
    colnames(dat_sub_area) <- c("subbas", "area")
    # merge data and calculate area-weighted mean and time period sum
    dat_surf_mm <- left_join(dat_tmp,
                           as.data.frame(dat_sub_area) %>%
                             mutate(subbas = as.integer(subbas), area_sum = sum(area)),
                           by = "subbas") %>%
      group_by(date, area_sum) %>%
      summarise(value_sum = sum(value)) %>% # sum over catchment (m3/timestep)
      mutate(value_mm = value_sum*1000 / (area_sum*1e6)) %>% # daily sum over catchment in (mm)
      ungroup() %>%
      summarise(value_sum_mm = sum(value_mm)) # sum over simulation period (mm)

    out[["runoff_surf_mm"]] <- dat_surf_mm$value_sum_mm

  }
  if("runoff_total_mm" %in% return_val) {
    dat_tmp <- filter(dat_wasa, group == "runoff_total_mm")
    # get catchment area (m2)
    dat_sub <- readLines(paste(dir_run, "input/Hillslope/hymo.dat", sep="/"))
    dat_sub <- dat_sub[-c(1,2)]
    dat_sub_area <- sapply(dat_sub, function(x) as.numeric(unlist(strsplit(x, "\t"))[c(1,2)]), USE.NAMES = F)
    dat_sub_area <- dat_sub_area[2, order(dat_sub_area[1,])]
    dat_sub_area <- sum(dat_sub_area) * 1e6
    # sum of catchment river outflow (mm)
    outflow <- dat_tmp$value * 60*60*24 # m3/day
    outflow <- sum(outflow) # m3
    outflow <- outflow *1000 # L
    outflow <- outflow / dat_sub_area # L/m2 = mm

    out[["runoff_total_mm"]] <- outflow

  }
  if(any(c("eta_mm", "etp_mm") %in% return_val)) {
    dat_tmp <- filter(dat_wasa, group %in% c("eta_mm", "etp_mm"))
    # get catchment area (m2)
    dat_sub <- readLines(paste(dir_run, "input/Hillslope/hymo.dat", sep="/"))
    dat_sub <- dat_sub[-c(1,2)]
    dat_sub_area <- t(sapply(dat_sub, function(x) as.numeric(unlist(strsplit(x, "\t"))[c(1,2)]), USE.NAMES = F))
    colnames(dat_sub_area) <- c("subbas", "area")
    # merge data and calculate area-weighted mean and time period sum
    dat_sums <- left_join(dat_tmp,
                         as.data.frame(dat_sub_area) %>%
                           mutate(subbas = as.integer(subbas), area_sum = sum(area), area_weight = area / area_sum),
                         by = "subbas") %>%
      mutate(value_w = value * area_weight) %>%
      group_by(group, date) %>%
      summarise(value = sum(value_w)) %>% # daily area-weighted means over catchment
      group_by(group) %>%
      summarise(value_sum = sum(value)) # sum over simulation period

    out[dat_sums$group] <- dat_sums$value_sum

  }

  if(length(out) > 1 && class(out) != "list") out <- as.list(out)

  # clean up
  if(!keep_rundir) unlink(dir_run, recursive = T)

  # output
  return(out)
} # EOF
