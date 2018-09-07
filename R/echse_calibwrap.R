#' Wrapper function for calibration of the WASA engine of the ECHSE environment
#'
#' A wrapper function, executing a simulation of the WASA engine for a given
#' parameter realisation.
#'
#' @param pars A named vector of type numeric with the values of selected parameters.
#' Replaces all specified parameter values within the sharedParamNum_WASA_* files
#' with the given values.
#'
#' @param dir_input Character string giving the main simulation directory containing a
#' directory 'data' with the readily prepared ECHSE input (e.g. prepared with function
#' \code{\link[WasaEchseTools]{echse_prep_runs}}).
#'
#' @param dir_run Character specifying the directory for the model run with the current
#' parameter realisation. Default: A temporary directory created with \code{\link{tempfile}}.
#'
#' @param echse_app Character string giving the system command of the application.
#'
#' @param choices A named data.frame with each element containing the flag for a specific
#' choice. See the latest version of ECHSE's WASA engine for required choice flags.
#' Replaces all specified parameter values within the sharedParamNum_WASA_* files
#' with the given values. If \code{NULL} (the default), it is assumed the model is
#' run with default or manually defined selections.
#'
#' @param sim_start Character string giving the start date of the simulation period
#' in the format "\%Y-\%m-\%d \%H:\%M:\%S".
#'
#' @param sim_end Character string giving the end date of the simulation period
#' in the format "\%Y-\%m-\%d \%H:\%M:\%S".
#'
#' @param resolution Integer giving the simulation time step in seconds. Default:
#' 86400 (i.e. daily resolution).
#'
#' @param dat_streamflow OPTIONAL: Object of class 'xts' containing a time series of
#' streamflow at the catchment outlet in (m3/s) in the resolution of the model run.
#' NA values are allowed and will be discarded for goodness of fit calculations.
#' Only needed if \code{return_val = 'nse'}.
#'
#' @param dat_pr OPTIONAL: Object of class 'xts' containing a time series of catchment-wide
#' average precipitation in (m3) in the resolution of the model run. If not given (default),
#' it will be read (and calculated) from the ECHSE input files if needed. However,
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
#' @param warmup_start Character string giving the start date of the warm-up period
#' in the format "\%Y-\%m-\%d \%H:\%M:\%S". If \code{NULL} (default), argument 'sim_start'
#' will be used.
#'
#' @param warmup_len Integer giving the length of the warm-up period in months. Default: 3.
#'
#' @param max_pre_runs Integer specifying the maximum number of warm-up iterations to be
#' applied. If the relative storage change is still larger than \code{storage_tolerance}
#' after \code{max_pre_runs} iterations, the warm-up will be aborted and the model be run
#' anyway. A warning will be issued. Default: 20.
#'
#' @param storage_tolerance Numeric value giving the relative change of the model's water
#' storages between two connsecutive warm-up runs below which the warm-up will be
#' concluded and the actual model simulation be started. Default: 0.01.
#'
#' @param return_val Character vector specifying your choice of what this function
#' shall return. Default: 'river_flow'. See description of return value below.
#'
#' @param keep_rundir Value of type \code{logical}. Shall directory \code{dir_run}
#' be retained (\code{TRUE}) or deleted (\code{FALSE}) after function execution?
#' Default: \code{FALSE}.
#'
#' @param error2warn Value of type \code{logical}. Shall runtime errors of the model be
#' reported as a warning instead of stopping this function with an error? In case
#' of reasonable model output, this wrapper function will proceed as usual. If no
#' reasonable output could be found, a value \code{NA} will be returned.
#' Default: \code{FALSE}.
#'
#' @param nthreads Number of cores that shall be employed for the ECHSE run (argument
#' 'number_of_threads' in configuration file). See ECHSE core manual for more information.
#'
#' @details The function can be employed by model calibration functions such as
#' \code{\link[HydroBayes]{dream}} or \code{\link[ppso]{optim_dds}}, or to execute
#' single model runs within a single function call.
#'
#' @note In the current form, this function cannot deal with *_tpl.dat files generated
#' by \code{\link[WasaEchseTools]{echse_prep_runs}} (argument \code{prep_tpl})!
#'
#' @return Function returns a named list or a single element (if \code{length(return_val) = 1})
#' of numeric value(s). Can be controlled by argument \code{return_val}. Currently
#' implemented are the options:
#'
#' river_flow: Simulated river flow leaving the catchment outlet in m3/s for the
#' specified simulation period and resolution
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
echse_calibwrap <- function(
  pars = NULL,
  dir_input = NULL,
  dir_run = paste0(tempfile(pattern = "echse_calib_"), sep="/"),
  echse_app = NULL,
  choices = NULL,
  sim_start = NULL,
  sim_end = NULL,
  resolution = 86400,
  dat_streamflow = NULL,
  dat_pr = NULL,
  flood_thresh = NULL,
  thresh_zero = NULL,
  warmup_start=NULL,
  warmup_len = 3,
  max_pre_runs = 20,
  storage_tolerance = 0.01,
  return_val = "river_flow",
  keep_rundir = FALSE,
  error2warn = FALSE,
  nthreads = 1
) {

  if("hydInd" %in% return_val) {
    if(!is.numeric(flood_thresh)) stop("Argument return_val = hydInd requires argument flood_thresh to be given!")
    if(!is.numeric(thresh_zero)) stop("Argument return_val = hydInd requires argument thresh_zero to be given!")
  }
  if("nse" %in% return_val) {
    if(!is.xts(dat_streamflow)) stop("Argument return_val = nse requires an object of class 'xts' for argument dat_streamflow!")
  }

  # create run directory
  run_pars <- paste(dir_run, "pars", sep="/")
  run_out <- paste(dir_run, "out", sep="/")
  dir.create(paste(dir_run, sep="/"), recursive = T, showWarnings = F)
  dir.create(run_pars, recursive = T, showWarnings = F)
  dir.create(run_out, recursive = T, showWarnings = F)

  # in some cases (if there is only one subbasin), there is no WASA_rch object -> special treatment needed
  dat_objdecl <- read.table(paste(dir_input, "data/catchment/objDecl.dat", sep="/"), sep="\t", header=T)
  if(!any(grepl("WASA_rch", dat_objdecl$objectGroup))) {
    rch_missing <- TRUE
  } else {
    rch_missing <- FALSE
  }

  sharedsvc_path <- paste(run_pars, "sharedParamNum_WASA_svc.dat", sep="/")
  sharedlu_path <- paste(run_pars, "sharedParamNum_WASA_lu.dat", sep="/")
  sharedrch_path <- paste(run_pars, "sharedParamNum_WASA_rch.dat", sep="/")

  # ADJUST PARAMETERS AND CHOICES #
  # replace log values
  if(any(grepl(names(pars), pattern = "^log_"))) {
    log_trans = grepl(names(pars), pattern = "^log_") #find log-transformed parameters
    pars[log_trans]=exp(pars[log_trans]) #transform back to non-log scale
    names(pars) = sub(names(pars), pattern = "^log_", rep="") #remove "log_" from name
  }
  # iterate over all sharedParamNum_*.dat files (in current implementation only parameters of type 'sharedParamNum' can be calibrated)
  par_files <- dir(paste(dir_input, "data/parameter/", sep="/"), pattern = "sharedParamNum_")
  #if(rch_missing) par_files <- subset(par_files, !grepl("WASA_rch", par_files))
  for(f in par_files) {
    # read parameter file and replace parameters
    dat <- read.table(paste(dir_input, "data/parameter", f, sep="/"), header = T, sep="\t", stringsAsFactors = FALSE) %>%
      left_join(.,
                data.frame(parameter = names(pars), val_repl = pars, stringsAsFactors = FALSE),
                by = "parameter") %>%
      mutate(value = ifelse(is.na(val_repl), value, val_repl)) %>%
      dplyr::select(-val_repl) %>%
      { # adjust choices if available
      if(!is.null(choices)) {
        left_join(.,
                  data.frame(choices = names(choices), val_repl = unlist(choices), stringsAsFactors = FALSE),
                  by = c("parameter" = "choices")) %>%
        mutate(value = ifelse(is.na(val_repl), value, val_repl)) %>%
        dplyr::select(-val_repl)
      } else .
      }

    # write output
    write.table(dat, paste(run_pars, f, sep="/"), quote=F, row.names=F, sep="\t")
  }

  # WARM-UP RUNS #
  # get subbasin parameters (needed to calculate storages)
  sub_pars <- read.table(paste(dir_input, "data", "parameter", "paramNum_WASA_sub.dat", sep="/"), header=T)

  # output debug file
  write(c("# Objects for which debug info is to be printed", "object", "none"),
        file = paste(run_pars, "output_debug.txt", sep="/"), sep="\n")

  # output selection file
  output_warmup <- data.frame(object = rep(sub_pars$object, each=2),
                              variable = rep(c("v_runstor", "v_soilwat"), times=nrow(sub_pars)),
                              digits = 12)
  write.table(output_warmup, paste(run_pars, "output_selection.txt", sep="/"), row.names = F, col.names = T, sep="\t", quote=F)

  # output state file
  write(c("time", "1970-01-01 00:00:00"), file = paste(run_pars, "output_state.txt", sep="/"), sep="\n")

  # copy model states into run par directory
  file.copy(paste(dir_input, "data/initials/init_scal.dat", sep="/"),paste(run_pars, "init_scal.dat", sep="/"), overwrite = T)
  file.copy(paste(dir_input, "data/initials/init_vect.dat", sep="/"),paste(run_pars, "init_vect.dat", sep="/"), overwrite = T)

  # adjust model config file
  if(is.null(warmup_start)) warmup_start <- sim_start
  warmup_end <- seq(as.POSIXct(warmup_start, tz='UTC'), by=paste(warmup_len, "month"), length=2)[2]-resolution
  warmup_end <- format(warmup_end, "%Y-%m-%d %H:%M:%S")
  model_cnf <- readLines(system.file("echse_ctrl_tpl/cnf_default", package="WasaEchseTools"))
  model_cnf <- gsub("NCORES",  nthreads, model_cnf)
  model_cnf <- gsub("MODELDIR", paste(dir_input, "data", sep="/"), model_cnf)
  model_cnf <- gsub("OUTDIR",  run_pars, model_cnf)
  model_cnf <- gsub("RUNSTART", warmup_start, model_cnf)
  model_cnf <- gsub("RUNEND", warmup_end, model_cnf)
  model_cnf <- gsub("RESOLUTION", resolution, model_cnf)
  model_cnf <- gsub("INITSCALFILE", paste(run_pars, "init_scal.dat", sep="/"), model_cnf)
  model_cnf <- gsub("INITVECTFILE", paste(run_pars, "init_vect.dat", sep="/"), model_cnf)
  model_cnf <- gsub("SHAREDPARSVC", sharedsvc_path, model_cnf)
  model_cnf <- gsub("SHAREDPARLU", sharedlu_path, model_cnf)
  model_cnf <- gsub("SHAREDPARRCH", sharedrch_path, model_cnf)
  if(rch_missing) {
    model_cnf <- gsub("paramNum_WASA_rch.dat", "dummy_num.dat", model_cnf)
    model_cnf <- gsub("paramFun_WASA_rch.dat", "dummy_fun.dat", model_cnf)
  }
  writeLines(model_cnf, paste(run_pars, "cnf_default", sep="/"))

  # model arguments
  model_args <- c(file_control=paste(run_pars, "cnf_default", sep="/"),
                  file_log=paste(run_out, "run.log", sep="/"),
                  file_err=paste(run_out, "run.err.html", sep="/"),
                  format_err="html",
                  silent="true",
                  outputDirectory=run_out)

  # construct system command to run model
  cmd <- paste(c(echse_app, paste(names(model_args), model_args, sep="=")), collapse = " ")

  # initialise water storage tracker
  storage_before <- 0

  # loop over pre-runs
  for (i in 1:max_pre_runs) {
    # run model
    status <- system(cmd, intern=FALSE, ignore.stderr=FALSE, wait=TRUE)
    if(file.exists(paste(run_out, "run.err.html", sep="/"))) {
      if(error2warn) {
        warning(paste0("ECHSE returned a runtime error during warm-up, see log file: ", paste(run_out, "run.err.html", sep="/"), ". Continue model run ..."))
        break
      } else {
        stop(paste("ECHSE returned a runtime error during warm-up, see log file:", paste(run_out, "run.err.html", sep="/")))
      }
    }

    # adjust state files
    file.copy(dir(run_out, "statesScal", full.names = T), paste(run_pars, "init_scal.dat", sep="/"), overwrite = T)
    file.copy(dir(run_out, "statesVect", full.names = T), paste(run_pars, "init_vect.dat", sep="/"), overwrite = T)

    # read in results and calculate current catchment-wide soil + groundwater storage
    storage_after <- sub_pars %>%
      # read subbasin output
      ddply("object", function(x) {
        read.table(paste0(run_out, "/", x$object, ".txt"), header=T, sep="\t") %>%
          dplyr::select(-end_of_interval) %>%
          filter(row_number()==n()) %>%
          mutate(area = x$area)
      }) %>%
      # tidy data.frame()
      melt(id.vars = c("object", "area"))  %>%
      # calculate catchment-wide water storage in (m3)
      mutate(storsum_t = value * area * 1e6) %>%
      summarise(storsum=sum(storsum_t))

    # clean output directory
    file.remove(dir(run_out, full.names = T))

    # compare with storage from previous run
    rel_storage_change <- abs(storage_after-storage_before)
    # avoid NaNs sum(storage_before)==0
    if (storage_before!=0) rel_storage_change <- rel_storage_change/storage_before

    # check if storage changes are below tolerance limit
    if (rel_storage_change < storage_tolerance) break
    storage_before <- storage_after
  }
  if(i == max_pre_runs)
    warning(paste("Relative storage change after 'max_pre_runs' iterations was still above the tolerance threshold: ", round(rel_storage_change, 3)))



  # ACTUAL MODEL RUN #
  # output debug file
  content <- c("# Objects for which debug info is to be printed", "object", "none")
  write(content, file = paste(run_pars, "output_debug.txt", sep="/"), sep="\n")

  # output selection file
  output_sel <- data.frame(object = "node_su_out_1", variable = "out", digits = 12)
  write.table(output_sel, paste(run_pars, "output_selection.txt", sep="/"), row.names = F, col.names = T, sep="\t", quote=F)

  # output state file
  write(c("time", "1970-01-01 00:00:00"), file = paste(run_pars, "output_state.txt", sep="/"), sep="\n")

  # adjust model config file
  model_cnf <- readLines(system.file("echse_ctrl_tpl/cnf_default", package="WasaEchseTools"))
  model_cnf <- gsub("NCORES",  nthreads, model_cnf)
  model_cnf <- gsub("MODELDIR", paste(dir_input, "data", sep="/"), model_cnf)
  model_cnf <- gsub("OUTDIR",  run_pars, model_cnf)
  model_cnf <- gsub("RUNSTART", sim_start, model_cnf)
  model_cnf <- gsub("RUNEND", sim_end, model_cnf)
  model_cnf <- gsub("RESOLUTION", resolution, model_cnf)
  model_cnf <- gsub("INITSCALFILE", paste(run_pars, "init_scal.dat", sep="/"), model_cnf)
  model_cnf <- gsub("INITVECTFILE", paste(run_pars, "init_vect.dat", sep="/"), model_cnf)
  model_cnf <- gsub("SHAREDPARSVC", sharedsvc_path, model_cnf)
  model_cnf <- gsub("SHAREDPARLU", sharedlu_path, model_cnf)
  model_cnf <- gsub("SHAREDPARRCH", sharedrch_path, model_cnf)
  if(rch_missing) {
    model_cnf <- gsub("paramNum_WASA_rch.dat", "dummy_num.dat", model_cnf)
    model_cnf <- gsub("paramFun_WASA_rch.dat", "dummy_fun.dat", model_cnf)
  }
  writeLines(model_cnf, paste(run_pars, "cnf_default", sep="/"))

  # model arguments
  model_args <- c(file_control=paste(run_pars, "cnf_default", sep="/"),
                  file_log=paste(run_out, "run.log", sep="/"),
                  file_err=paste(run_out, "run.err.html", sep="/"),
                  format_err="html",
                  silent="true",
                  outputDirectory=run_out)

  # construct system command to run model
  cmd <- paste(c(echse_app, paste(names(model_args), model_args, sep="=")), collapse = " ")

  # run model
  status <- system(cmd, intern=FALSE, ignore.stderr=FALSE, wait=TRUE)
  if(file.exists(paste(run_out, "run.err.html", sep="/"))) {
    if(error2warn) {
      warning(paste0("ECHSE returned a runtime error during simulation, see log file: ",
                     paste(run_out, "run.err.html", sep="/"), ". Keep on running, returning 'NA'."))
      return(NA)
    } else {
      stop(paste("ECHSE returned a runtime error during simulation, see log file:", paste(run_out, "run.err.html", sep="/")))
    }
  }

  # get simulations
  file_echse <- paste(run_out, "node_su_out_1.txt", sep="/")
  tryCatch(dat_echse <- read.table(file_echse, header=T, sep="\t"),
           error = function(e) stop(paste("Error reading", file_echse, ":", e)))
  dat_echse <- dat_echse %>%
    mutate(date = as.POSIXct(.[[1]], tz ="UTC")-resolution, group = "echse") %>% # convert date to "begin of interval" (as in WASA output)
    rename(value = "out") %>%
    dplyr::select(date, group, value)

  # ignore errors if desired
  if(nrow(dat_echse) == 0) {
    if(!error2warn) {
      stop(paste0("There could be no output extracted from a model run, see ", run_out))
    } else {
      warning(paste0("No reasonable model output could be extracted from model run in ", run_out, ". NA will be returned!"))
      return(NA)
    }
  }

  # prepare output according to specifications
  out <- NULL
  if("hydInd" %in% return_val) {
    dat_sim_xts <- xts(dat_echse$value, dat_echse$date)
    # precipitation (model forcing); get catchment-wide value (area-weighted precipitation mean)
    if(is.null(dat_pr)) {
      datafiles <- read.table(paste(dir_input, "data/forcing/inputs_ext_datafiles.dat", sep="/"), sep="\t", header=T)
      prec_file <- as.character(datafiles$file[grep("precip", datafiles$variable)])
      prec_dat <- read.table(prec_file, header=T, sep="\t")
      dat_wgt <- read.table(paste(dir_input, "data/forcing/inputs_ext_locations.dat", sep="/"), sep="\t", header=T)
      dat_sub <- read.table(paste(dir_input, "data/parameter/paramNum_WASA_sub.dat", sep="/"), header=T)
      dat_prec <- left_join(dat_wgt %>%
                         filter(variable == "precip") %>%
                         dplyr::select(-variable) %>%
                         mutate_if(is.factor, as.character) %>%
                         separate(object, c("dum1", "sub"), sep="_", extra = "drop") %>%
                         dplyr::select(-dum1) %>%
                         distinct(),
                       prec_dat %>%
                         mutate(datetime = as.POSIXct(datetime, tz="UTC")) %>%
                         filter(datetime >= as.POSIXct(sim_start, tz ="UTC") & datetime <= as.POSIXct(sim_end, tz ="UTC")) %>%
                         melt(id.var = "datetime", variable.name = "location", value.name = "value") %>%
                         mutate_if(is.factor, as.character),
                       by = "location") %>%
        # calculate precipitation for every subbasin; in case of hourly resolution, calculate daily sums (hydInd() expects daily values)
        mutate(datetime = format(datetime, "%Y-%m-%d")) %>%
        group_by(sub, location, weight, datetime) %>%
        summarise(value = sum(value)) %>%
        group_by(sub, datetime) %>%
        summarise(value = sum(value*weight)) %>%
        ungroup() %>%
        # calculate catchment-wide precipitation in m3
        left_join(x = mutate(., sub = paste("sub", sub, sep="_")),
                  y = dat_sub %>% mutate_if(is.factor, as.character),
                  by = c("sub" = "object")) %>%
        group_by(datetime, sub, area) %>%
        summarise(value = sum(value)) %>%
        group_by(datetime) %>%
        # sums over the catchment in m3
        summarise(value = sum(value * area*1e3))
      # xts object
      dat_pr <- xts(dat_prec$value, as.POSIXct(dat_prec$datetime, tz ="UTC"))
    }
    # calculate diagnostic values
    out_vals <- suppressWarnings(hydInd(dat_sim_xts, dat_pr, na.rm = T, thresh.zero = thresh_zero, flood.thresh = flood_thresh))
    out_vals[which(is.na(out_vals) | is.nan(out_vals))] <- 0

    out[["hydInd"]] <- out_vals

  }
  if("river_flow" %in% return_val) {
    out_vals <- xts(dat_echse$value, dat_echse$date)

    out[["river_flow"]] <- out_vals

  }
  if("nse" %in% return_val) {
    dat_nse <- left_join(select(dat_echse, date, value),
                         data.frame(obs = dat_streamflow,
                                    date = index(dat_streamflow)),
                         by = "date") %>%
      mutate(diffsq = (value - obs)^2,
             obsmean = mean(obs), diffsqobs = (obs - obsmean)^2) %>%
      summarise(nse = 1 - sum(diffsq) / sum(diffsqobs))
    out_vals <- as.numeric(dat_nse)

    out[["nse"]] <- out_vals

  }

  # clean up
  if(!keep_rundir) unlink(paste(dir_run), recursive = T)

  # output
  return(out)

} # EOF
