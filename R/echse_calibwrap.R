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
#' @param keep_warmup_states Logical value. Shall state file after finishing the warm-up
#' runs be stored for upcomming runs? In addition, file 'init_stor_sum_m3.dat' containg
#' a single value of total catchment storage in m3 will be created. WARNING: The files
#' will be stored in the ECHSE setup in \code{dir_input} and existing state files will be overwritten!
#' Default: \code{FALSE}. This option might be useful for calibration runs as it
#' might reduce the number of necessary warm-up runs for each new parameter set.
#'
#' @param return_val Character vector specifying your choice of what this function
#' shall return. Default: 'river_flow'. See description of return value below.
#'
#' @param return_sp Logical flag specifying if certain output shall be given including
#' spatial variability at subbasin level instead of mere catchment specific values.
#' See description of return values below for spatial outputs. Default: \code{FALSE}.
#'
#' @param log Character containg the name (and path) of a file into which information
#' about the function call will be written. See below for more information. Default:
#' \code{NULL}, i.e. no logile will be created. If the file already exists, the file
#' name to be created will be extended by a call to \code{\link{tempfile}}.
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
#' @return Function returns a named list or a single element (if \code{length(return_val) = 1})
#' of numeric value(s). Can be controlled by argument \code{return_val}. Currently
#' implemented are the options:
#'
#' river_flow_[mm|ts]: A single value of the catchment's simulated river outflow over
#' the specified simulation period in (mm) or a time series of class 'xts' in (m3/s).
#'
#' eta_mm[_ts]: A single value of the catchment's simulated amount of actual evapotranspiration
#' over the specified simulation period in (mm) or a time series of class 'xts' in (mm/timestep).
#' Respects \code{return_sp}.
#'
#' etp_mm[_ts]: A single value of the catchment's simulated amount of potential evapotranspiration
#' over the specified simulation period in (mm) or a time series of class 'xts' in (mm/timestep).
#' Respects \code{return_sp}.
#'
#' runoff_total_mm[_ts]: A single value of the catchment's simulated amount of total runoff,
#' i.e. runoff contribution into river(s), over the specified simulation period in (mm).
#' Respects \code{return_sp} or a time series of class 'xts' in (mm/timestep).
#'
#' runoff_gw_mm[_ts]: A single value of the catchment's simulated amount of groundwater runoff,
#' i.e. groundwater contribution into river(s), over the specified simulation period in (mm).
#' Respects \code{return_sp} or a time series of class 'xts' in (mm/timestep).
#'
#' runoff_sub_mm[_ts]: A single value of the catchment's simulated amount of subsurface runoff,
#' i.e. near-surface soil water (excl. groundwater!) contribution into river(s), over the specified
#' simulation period in (mm) or a time series of class 'xts' in (mm/timestep). Note: due to computational reasons, 'runoff_gw_mm' will
#' be given in addition if this value is set.
#' Respects \code{return_sp}.
#'
#' runoff_surf_mm[_ts]: A single value of the catchment's simulated amount of surface runoff,
#' i.e. surface water contribution into river(s), over the specified simulation period in (mm)
#' or a time series of class 'xts' in (mm/timestep).
#' Respects \code{return_sp}.
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
#'
#' If argument \code{log} was given, the logfile contains the following information:
#' parameter names and corresponding realisations, output element names and corresponding
#' values, \code{run_dir}: realisation of argument \code{dir_run}, \code{time_total}:
#' total time for function execution (i.e. until writing logfile, in secs.), \code{time_simrun}:
#' runtime of the simulation run (in secs.), \code{time_warmup}: total runtime for all
#' warm-up runs (in secs.), \code{warmup_iterations}: number of warm-up iterations,
#' \code{warmup_storchange}: relative storage change after the last warm-up iteration.
#'
#'
#' @note In the current form, this function cannot deal with *_tpl.dat files generated
#' by \code{\link[WasaEchseTools]{echse_prep_runs}} (argument \code{prep_tpl})!
#'
#' To avoid warm-up runs, set \code{max_pre_runs} or \code{warmup_len} to zero.
#'
#' \strong{Model's water balance}: runoff_total_mm = runoff_surf_mm + runoff_sub_mm + runoff_gw_mm.
#' Moreover, river_flow_mm should be slight less than (due to transmission losses and storage
#' changes; but all in all almost equal to) runoff_total_mm as long as there is no
#' artificial influence (e.g. from reservoirs).
#'
#' \strong{Initial states} can have a large influence, depending on setup and parameterisation.
#' Therefore, it is advisable to make some test runs with different warmup parameters (
#' \code{warmup_start}, \code{warmup_len}, \code{max_pre_runs}, \code{storage_tolerance},
#' \code{keep_warmup_states}) before starting a calibration run! Consecutive runs with
#' the same parameterisation should come up with the same results when \code{keep_warmup_states = TRUE}.
#'
#' If running this function in \strong{parallel} mode:
#' Setting argument \code{keep_warmup_states = TRUE} can have problematic side effects
#' due to parallel reading and (over-)writing of the same files. To accelerate
#' warmup anyway, you can generate default initial storage files by running a
#' default parametrisation over the warmup horizon until convergence is reached
#' and use the resulting storage files as initial storages for all calibration runs.
#' A similar problem arises When specifying argument \code{log}. In that case,
#' create an empty initial logfile '\code{log}' before starting the parallel
#' calibration routine to invoke the file name extension via \code{\link{tempfile}}
#' right away.
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
  keep_warmup_states = FALSE,
  return_val = "river_flow",
  return_sp = FALSE,
  log = NULL,
  keep_rundir = FALSE,
  error2warn = FALSE,
  nthreads = 1
) {
  time_start <- Sys.time()
  return_na <- FALSE # set to true in tryCatch(...) if error2warn = T
  if("hydInd" %in% return_val) {
    if(!is.numeric(flood_thresh)) stop("Argument return_val = hydInd requires argument flood_thresh to be given!")
    if(!is.numeric(thresh_zero)) stop("Argument return_val = hydInd requires argument thresh_zero to be given!")
  }
  if(any(c("nse", "kge") %in% return_val)) {
    if(!is.xts(dat_streamflow)) stop("Argument return_val = nse requires an object of class 'xts' for argument dat_streamflow!")
  }
  if(is.null(sim_start) | is.null(sim_end)) stop("Arguments 'sim_start' and 'sim_end' must be given!")
  tryCatch({
    sim_start <- as.POSIXct(sim_start, tz='UTC')
    sim_end <- as.POSIXct(sim_end, tz='UTC') + resolution # +1 time step due to date conversions in ECHSE ("bing of interval" vs. "end of interval" issue)
  }, error = function(e) stop("A problem occurred when coercing 'sim_start' and/or 'sim_end' to POSIX oject!"))
  if(!is.null(log)) {
    f_log <- TRUE
    logfile = log
    logdir <- sub(basename(logfile), "", logfile)
    if(file.exists(logfile)) logfile <- tempfile(sub(".[a-z]+$", "", basename(logfile)), tmpdir = logdir, fileext = sub("^[a-zA-Z0-9_-]+", "", basename(logfile)))
    if(logdir != "") dir.create(logdir, recursive = T, showWarnings = F)
    tryCatch(write(NULL, logfile), error = function(e) stop(paste("Could not create the given log file:", e)))
  } else {
    f_log <- FALSE
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

  #### ADJUST PARAMETERS AND CHOICES ####
  # replace log values
  if(any(grepl(names(pars), pattern = "^log_"))) {
    log_trans = grepl(names(pars), pattern = "^log_") #find log-transformed parameters
    pars_t <- pars # keep param names and values incl. "log_" for logfiles
    pars_t[log_trans]=exp(pars_t[log_trans]) #transform back to non-log scale
    names(pars_t) = sub(names(pars_t), pattern = "^log_", rep="") #remove "log_" from name
  }
  # iterate over all sharedParamNum_*.dat files (in current implementation only parameters of type 'sharedParamNum' can be calibrated)
  par_files <- dir(paste(dir_input, "data/parameter/", sep="/"), pattern = "sharedParamNum_")
  #if(rch_missing) par_files <- subset(par_files, !grepl("WASA_rch", par_files))
  for(f in par_files) {
    # read parameter file and replace parameters
    dat <- read.table(paste(dir_input, "data/parameter", f, sep="/"), header = T, sep="\t", stringsAsFactors = FALSE) %>%
      left_join(.,
                data.frame(parameter = names(pars_t), val_repl = pars_t, stringsAsFactors = FALSE),
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

  # get subbasin parameters (needed to calculate storages and certain outputs)
  sub_pars <- read.table(paste(dir_input, "data", "parameter", "paramNum_WASA_sub.dat", sep="/"), header=T, stringsAsFactors = F)


  #### WARM-UP RUNS ####

  # copy model states into run par directory
  file.copy(paste(dir_input, "data/initials/init_scal.dat", sep="/"),paste(run_pars, "init_scal.dat", sep="/"), overwrite = T)
  file.copy(paste(dir_input, "data/initials/init_vect.dat", sep="/"),paste(run_pars, "init_vect.dat", sep="/"), overwrite = T)

  # conduct warm-up runs?
  time_warmup <- NA
  i_warmup <- NA
  rel_storage_change <- NA
  if(warmup_len > 0 & max_pre_runs > 0) {
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

    # adjust model config file
    if(is.null(warmup_start)) warmup_start <- sim_start
    warmup_start <- as.POSIXct(warmup_start, tz="UTC")
    warmup_end <- seq(warmup_start, by=paste(warmup_len, "month"), length=2)[2]-resolution
    warmup_end <- format(warmup_end, "%Y-%m-%d %H:%M:%S")
    model_cnf <- readLines(system.file("echse_ctrl_tpl/cnf_default", package="WasaEchseTools"))
    model_cnf <- gsub("NCORES",  nthreads, model_cnf)
    model_cnf <- gsub("MODELDIR", paste(dir_input, "data", sep="/"), model_cnf)
    model_cnf <- gsub("OUTDIR",  run_pars, model_cnf)
    model_cnf <- gsub("RUNSTART", format(warmup_start, "%Y-%m-%d %H:%M:%S"), model_cnf)
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
    if(file.exists(paste(dir_input, "data/initials/init_stor_sum_m3.dat", sep="/"))) {
      storage_before <- as.numeric(read_file(paste(dir_input, "data/initials/init_stor_sum_m3.dat", sep="/")))
    } else {
      storage_before <- 0
    }

    # loop over pre-runs
    time_warmup <- system.time({
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
          bind_cols(.$object %>%
                      map_dfr(function(x) {
                                suppressMessages(read_tsv(paste0(run_out, "/", x, ".txt"))) %>%
                                  select(-end_of_interval) %>%
                                  filter(row_number()==n())
                              }),
                    .) %>%
          gather(key = variab, value = value, -object, -area) %>%
          # calculate catchment-wide water storage in (m3)
          mutate(storsum_t = value * area * 1e6) %>%
          summarise(storsum=sum(storsum_t))
        storage_after <- storage_after$storsum

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
    }) # measure warmup time
    i_warmup <- i
    if(i_warmup == max_pre_runs)
      warning(paste("Relative storage change after 'max_pre_runs' iterations was still above the tolerance threshold: ", round(rel_storage_change, 3)))

    # replace old state file in dir_input?
    if(keep_warmup_states) {
      file.copy(paste(run_pars, "init_scal.dat", sep="/"), paste(dir_input, "data/initials/init_scal.dat", sep="/"), overwrite = T)
      file.copy(paste(run_pars, "init_vect.dat", sep="/"), paste(dir_input, "data/initials/init_vect.dat", sep="/"), overwrite = T)
      write(storage_after, file=paste(dir_input, "data/initials/init_stor_sum_m3.dat", sep="/")) # tedious to calculate storage sum from init_*.dat files
    }

  } # warmup run to be conducted?


  #### ACTUAL MODEL RUN ####
  # output debug file
  content <- c("# Objects for which debug info is to be printed", "object", "none")
  write(content, file = paste(run_pars, "output_debug.txt", sep="/"), sep="\n")

  # output selection file
  output_sel <- data.frame(object = NULL, variable = NULL, digits = NULL)
  if(any(c("nse", "kge", "hydInd", "river_flow_ts", "river_flow_mm") %in% return_val))
    output_sel <- bind_rows(output_sel, data.frame(object = "node_su_out_1", variable = "out", digits = 12, stringsAsFactors = F))

  if(any(str_detect(return_val, "runoff_surf")))
    output_sel <- bind_rows(output_sel, data.frame(object = sub_pars$object, variable = "r_out_surf", digits = 12, stringsAsFactors = F))
  if(any(str_detect(return_val, "runoff_sub")))
    output_sel <- bind_rows(output_sel, data.frame(object = sub_pars$object, variable = "r_out_inter", digits = 12, stringsAsFactors = F))
  if(any(str_detect(return_val, "runoff_gw")))
    output_sel <- bind_rows(output_sel, data.frame(object = sub_pars$object, variable = "r_out_base", digits = 12, stringsAsFactors = F))
  if(any(str_detect(return_val, "runoff_total")))
    output_sel <- bind_rows(output_sel, data.frame(object = sub_pars$object, variable = "r_out_total", digits = 12, stringsAsFactors = F))
  return_grp <- c("runoff_total" = "r_out_total",
                  "runoff_surf" = "r_out_surf",
                  "runoff_sub" = "r_out_inter",
                  "runoff_gw" = "r_out_base")

  if(any(str_detect(return_val, "eta")))
    output_sel <- bind_rows(output_sel, data.frame(object = rep(sub_pars$object,2),
                                                   variable = rep(c("eta", "eti"), each=length(sub_pars$object)),
                                                   digits = 12, stringsAsFactors = F))
  if(any(str_detect(return_val, "etp")))
    output_sel <- bind_rows(output_sel, data.frame(object = sub_pars$object, variable = "etp", digits = 12, stringsAsFactors = F))

  write.table(output_sel, paste(run_pars, "output_selection.txt", sep="/"), row.names = F, col.names = T, sep="\t", quote=F)

  # output state file
  write(c("time", "1970-01-01 00:00:00"), file = paste(run_pars, "output_state.txt", sep="/"), sep="\n")

  # adjust model config file
  model_cnf <- readLines(system.file("echse_ctrl_tpl/cnf_default", package="WasaEchseTools"))
  model_cnf <- gsub("NCORES",  nthreads, model_cnf)
  model_cnf <- gsub("MODELDIR", paste(dir_input, "data", sep="/"), model_cnf)
  model_cnf <- gsub("OUTDIR",  run_pars, model_cnf)
  model_cnf <- gsub("RUNSTART", format(sim_start, "%Y-%m-%d %H:%M:%S"), model_cnf)
  model_cnf <- gsub("RUNEND", format(sim_end, "%Y-%m-%d %H:%M:%S"), model_cnf)
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
  time_simrun <- system.time({
    status <- system(cmd, intern=FALSE, ignore.stderr=FALSE, wait=TRUE)
  })
  if(file.exists(paste(run_out, "run.err.html", sep="/"))) {
    # write logfile
    if(f_log) {
      # meta information of model call
      time_end <- Sys.time()
      out_log <- data.frame(group = c(rep("pars", length(pars)), rep("meta", 6)),
                            variable = c(names(pars), "run_dir", "time_total", "time_simrun", "time_warmup", "warmup_iterations", "warmup_storchange"),
                            value = c(round(pars, 4), dir_run, round(difftime(time_end, time_start, units = "s"), 1), round(time_simrun["elapsed"], 1), round(time_warmup["elapsed"], 1), i_warmup, round(rel_storage_change, 5))
      )
      write.table(out_log, file=paste0(logfile, ".err"), sep="\t", quote=F, row.names=F, col.names=T)
    }
    if(error2warn) {
      warning(paste0("ECHSE returned a runtime error during simulation, see log file: ",
                     paste(run_out, "run.err.html", sep="/"), ". Keep on running, returning 'NA'."))
      return(NA)
    } else {
      stop(paste("ECHSE returned a runtime error during simulation, see log file:", paste(run_out, "run.err.html", sep="/")))
    }
  }

  # get simulations
  dat_echse <- NULL
  if("node_su_out_1" %in% output_sel$object) {
    file_echse <- paste(run_out, "node_su_out_1.txt", sep="/")
    tryCatch(dat_tmp <- read.table(file_echse, header=T, sep="\t"),
             error = function(e) stop(paste("Error reading", file_echse, ":", e)))
    dat_echse <- dat_tmp %>%
      mutate(date = as.POSIXct(.[[1]], tz ="UTC")-resolution, group = "river_flow", subbas = "all") %>% # convert date to "begin of interval" (as in WASA output)
      rename(value = "out") %>%
      dplyr::select(date, group, value) %>%
      bind_rows(dat_echse)
  }
  if(any(grepl("sub_[0-9]+", output_sel$object))) {
    dat_echse <- dir(run_out, "sub_[0-9]+", full.names = T) %>%
      map_dfr(function(f) {
        suppressMessages(read_tsv(f)) %>%
          mutate(file = gsub(".txt$", "", basename(f)))
      }) %>%
      mutate(date = as.POSIXct(.[[1]], tz ="UTC")-resolution) %>% # convert date to "begin of interval" (as in WASA output)
      select(-end_of_interval) %>%
      gather(key = group, value = value, -file, - date) %>%
      bind_rows(dat_echse)
  }

  # ignore errors if desired
  if(nrow(dat_echse) == 0) {
    # write logfile
    if(f_log) {
      # meta information of model call
      time_end <- Sys.time()
      out_log <- data.frame(group = c(rep("pars", length(pars)), rep("meta", 6)),
                            variable = c(names(pars), "run_dir", "time_total", "time_simrun", "time_warmup", "warmup_iterations", "warmup_storchange"),
                            value = c(round(pars, 4), dir_run, round(difftime(time_end, time_start, units = "s"), 1), round(time_simrun["elapsed"], 1), round(time_warmup["elapsed"], 1), i_warmup, round(rel_storage_change, 5))
      )
      write.table(out_log, file=paste0(logfile, ".err"), sep="\t", quote=F, row.names=F, col.names=T)
    }
    if(!error2warn) {
      stop(paste0("There could be no output extracted from a model run, see ", run_out))
    } else {
      warning(paste0("No reasonable model output could be extracted from model run in ", run_out, ". NA will be returned!"))
      return(NA)
    }
  }

  # prepare output according to specifications
  out <- NULL
  tryCatch({
    if(any(str_detect(return_val, "hydInd"))) {
      dat_tmp <- filter(dat_echse, group == "river_flow")
      dat_sim_xts <- xts(dat_tmp$value, dat_tmp$date)
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
                           filter(datetime >= sim_start & datetime <= sim_end) %>%
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

      out[names(out_vals)] <- out_vals
      out <- as.list(out)

    }
    if(any(str_detect(return_val, "river_flow"))) {
      out_vals_t <- dat_echse %>%
        filter(group == "river_flow")

      if("river_flow_ts" %in% return_val) {
        out_vals_ts <- xts(out_vals_t$value, out_vals_t$date)
        colnames(out_vals_ts) <- "river_flow"
        if(any(names(out) == "xts")) {
          out[["xts"]] <- cbind(out$xts, out_vals_ts)
        } else {
          out[["xts"]] <- out_vals_ts
        }
      }

      if("river_flow_mm" %in% return_val) {
        out_vals_mm <- out_vals_t %>%
          mutate(value = value * resolution) %>% # m3/timestep
          summarise(value = sum(value) * 1000 / (sum(sub_pars$area)*1e6))
        out[["river_flow_mm"]] <- out_vals$value
      }
    }
    if("nse" %in% return_val) {
      dat_tmp <- filter(dat_echse, group == "river_flow")
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
      dat_tmp <- filter(dat_echse, group == "river_flow")
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
    if(any(str_detect(return_val, "eta"))) {
      out_vals <- dat_echse %>%
        filter(group %in% c("eta", "eti")) %>%
        left_join(.,
                  sub_pars %>%
                    mutate(area_sum = sum(area), wgt = area/area_sum),
                  by = c("file" = "object")) %>%
        mutate(value = value * resolution*1000) %>% # unit mm/timestep
        group_by(file, wgt, date) %>%
        summarise(value = sum(value)) # sums of eta (i.e. eta total = eta + eti) over simulation period (mm)

      if("eta_mm" %in% return_val) {
        out_vals_agg <- out_vals %>%
          group_by(file, wgt) %>%
          summarise(value = sum(value)) %>%
          ungroup() %>%
          mutate(value_catch = sum(value * wgt))
        out[["eta_mm"]] <- unique(out_vals_agg$value_catch)
        if(return_sp) out[paste("eta_mm", out_vals_agg$file, sep="_")] <- out_vals_agg$value
      }

      if("eta_mm_ts" %in% return_val) {
        out_vals_catch <- out_vals %>%
          group_by(date) %>%
          summarise(value = sum(value * wgt)) # area-weighted sums of eta per timestep over whole catchment
        out_vals_catch <- xts(out_vals_catch$value, out_vals_catch$date)
        colnames(out_vals_catch) <- "eta_mm"
        if(any(names(out) == "xts")) {
          out[["xts"]] <- cbind(out$xts, out_vals_catch)
        } else {
          out[["xts"]] <- out_vals_catch
        }

        if(return_sp) {
          out_vals_sp <- out_vals %>%
            ungroup() %>%
            select(-wgt) %>%
            spread(key = file, value = value) %>%
            rename_at(vars(-date), funs(paste0("eta_mm_", .)))
          out_vals_sp <- xts(select(out_vals_sp, -date), out_vals_sp$date)
          out[["xts"]] <- cbind(out$xts, out_vals_sp)
        }
      }
    }
    if(any(str_detect(return_val, "etp"))) {
      out_vals <- dat_echse %>%
        filter(group == "etp") %>%
        left_join(.,
                  sub_pars %>%
                    mutate(area_sum = sum(area), wgt = area/area_sum),
                  by = c("file" = "object")) %>%
        mutate(value = value * resolution*1000) %>% # unit mm/timestep
        group_by(file, wgt, date)

      if("etp_mm" %in% return_val) {
        out_vals_agg <- out_vals %>%
          group_by(file, wgt) %>%
          summarise(value = sum(value)) %>%
          ungroup() %>%
          mutate(value_catch = sum(value * wgt))
        out[["etp_mm"]] <- unique(out_vals_agg$value_catch)
        if(return_sp) out[paste("etp_mm", out_vals_agg$file, sep="_")] <- out_vals_agg$value
      }

      if("etp_mm_ts" %in% return_val) {
        out_vals_catch <- out_vals %>%
          group_by(date) %>%
          summarise(value = sum(value * wgt)) # area-weighted sums of eta per timestep over whole catchment
        out_vals_catch <- xts(out_vals_catch$value, out_vals_catch$date)
        colnames(out_vals_catch) <- "etp_mm"
        if(any(names(out) == "xts")) {
          out[["xts"]] <- cbind(out$xts, out_vals_catch)
        } else {
          out[["xts"]] <- out_vals_catch
        }

        if(return_sp) {
          out_vals_sp <- out_vals %>%
            ungroup() %>%
            select(-wgt, -area, -area_sum, -group) %>%
            spread(key = file, value = value) %>%
            rename_at(vars(-date), funs(paste0("etp_mm_", .)))
          out_vals_sp <- xts(select(out_vals_sp, -date), out_vals_sp$date)
          out[["xts"]] <- cbind(out$xts, out_vals_sp)
        }
      }
    }
    if(any(str_detect(return_val, paste(names(return_grp), collapse = "|")))) {
      out_vals <- dat_echse %>%
        filter(group %in% return_grp) %>%
        left_join(.,
                  sub_pars %>%
                    mutate(area_sum = sum(area), wgt = area/area_sum),
                  by = c("file" = "object")) %>%
        mutate(value = value * resolution, # unit m3/timestep
               group = factor(group, levels = return_grp, labels = names(return_grp))) %>%
        mutate(value = value / (area*1000)) # unit mm/timestep

      # return aggregated values
      return_agg <- return_val[which(return_val %in% paste(names(return_grp), "mm", sep="_"))]
      if(length(return_agg) > 0) {
        out_vals_agg <- out_vals %>%
          filter(group %in% str_replace_all(return_agg, "_mm", "")) %>%
          group_by(group, file, wgt) %>%
          summarise(value = sum(value)) %>% # mm/timestep per subbasin for each variable
          group_by(group) %>%
          mutate(value_catch = sum(value*wgt)) # area-weightes sum over catchment (mm/timestep)
        out[paste(unique(out_vals_agg$group), "mm", sep="_")] <- unique(out_vals_agg$value_catch)
        if(return_sp)  out[paste(out_vals_agg$group, out_vals_agg$file, sep="_")] <- out_vals_agg$value
      }

      return_ts <- return_val[which(return_val %in% paste(names(return_grp), "mm_ts", sep="_"))]
      if(length(return_ts) > 0) {
        out_vals_catch <- out_vals %>%
          filter(group %in% str_replace_all(return_ts, "_mm_ts", "")) %>%
          group_by(group, date) %>%
          summarise(value = sum(value*wgt)) %>% # area-weighted catchment values in mm/timestep
          spread(key = group, value = value) %>%
          rename_at(vars(-date), funs(paste0(., "_mm")))
        out_vals_catch <- xts(select(out_vals_catch, -date), out_vals_catch$date)
        if(any(names(out) == "xts")) {
          out[["xts"]] <- cbind(out$xts, out_vals_catch)
        } else {
          out[["xts"]] <- out_vals_catch
        }

        if(return_sp) {
          out_vals_sp <- out_vals %>%
            filter(group %in% str_replace_all(return_ts, "_mm_ts", "")) %>%
            select(-wgt, -area, -area_sum) %>%
            unite(col="var", group, file, sep="_mm_") %>%
            spread(key = var, value = value)
          out_vals_sp <- xts(select(out_vals_sp, -date), out_vals_sp$date)
          out[["xts"]] <- cbind(out$xts, out_vals_sp)
        }
      }
    }

    if(length(out) > 1 && class(out) != "list") out <- as.list(out)

  }, error = function(e) {
    # write logfile
    if(f_log) {
      # meta information of model call
      time_end <- Sys.time()
      out_log <- data.frame(group = c(rep("pars", length(pars)), rep("meta", 6)),
                            variable = c(names(pars), "run_dir", "time_total", "time_simrun", "time_warmup", "warmup_iterations", "warmup_storchange"),
                            value = c(round(pars, 4), dir_run, round(difftime(time_end, time_start, units = "s"), 1), round(time_simrun["elapsed"], 1), round(time_warmup["elapsed"], 1), i_warmup, round(rel_storage_change, 5))
      )
      write.table(out_log, file=paste0(logfile, ".err"), sep="\t", quote=F, row.names=F, col.names=T)
    }
    if(!error2warn) {
      stop(paste0("Error while compiling output for model run, dir_run = ", dir_run))
    } else {
      warning(paste0("Could not compile output for model run in ", dir_run, ". NA will be returned!"))
      return_na <- TRUE
    }
  })
  if(return_na) return(NA)

  # clean up
  if(!keep_rundir) unlink(paste(dir_run), recursive = T)

  # write logfile
  if(f_log) {
    time_end <- Sys.time()
    if(any(names(out) == "xts")) {
      logfile_xts <- str_replace(logfile, ".dat$", "_xts.dat")
      write.zoo(out$xts, file=logfile_xts, index.name = "date", row.names = F, col.names = T, sep="\t", quote=F)
      outdat_log <- out[-grep("xts", names(out))]
    }
    if(length(outdat_log) > 0) {
      out_log <- data.frame(group = c(rep("pars", length(pars)), rep("output", length(outdat_log)), rep("meta", 6)),
                            variable = c(names(pars), names(outdat_log), "run_dir", "time_total", "time_simrun", "time_warmup", "warmup_iterations", "warmup_storchange"),
                            value = c(round(pars, 4), round(unlist(outdat_log),4), dir_run, round(difftime(time_end, time_start, units = "s"), 1), round(time_simrun["elapsed"], 1), round(time_warmup["elapsed"], 1), i_warmup, round(rel_storage_change, 5))
      )
    } else {
      out_log <- data.frame(group = c(rep("pars", length(pars)), rep("meta", 6)),
                            variable = c(names(pars), "run_dir", "time_total", "time_simrun", "time_warmup", "warmup_iterations", "warmup_storchange"),
                            value = c(round(pars, 4), dir_run, round(difftime(time_end, time_start, units = "s"), 1), round(time_simrun["elapsed"], 1), round(time_warmup["elapsed"], 1), i_warmup, round(rel_storage_change, 5))
      )
    }
    write.table(out_log, file=logfile, sep="\t", quote=F, row.names=F, col.names=T)
  }

  # output
  return(out)

} # EOF
