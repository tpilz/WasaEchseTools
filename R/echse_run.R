#' Executes simulations with the WASA engine of the ECHSE simulation environment
#'
#' Function executes simulation with the WASA engine of the ECHSE simulation environment
#' including warm-up runs to start from stable model states as initial conditions.
#'
#' @param dir_input Character string giving the main simulation directory containing a
#' directory 'data' with the readily prepared ECHSE input (e.g. prepared with function
#' \code{\link[WasaEchseTools]{echse_prep_runs}}).
#'
#' @param run_name Character string defining the name of a sub-directory within \code{dir_input}
#' to be created as storage of the actual simulations. It will further contain the directories
#' 'pars' containing ECHSE's control files and parameter files with the model structure choices
#' (if \code{choices} is given), and 'out' containing the actual simulation output.
#'
#' @param echse_app Character string giving the system command of the application.
#'
#' @param choices A named data.frame with each element containing the flag for a specific
#' choice. See the latest version of ECHSE's WASA engine for required choice flags. If
#' \code{NULL} (the default), it is assumed the model is run with default or manually
#' defined selections. In that case, the 'data/parameter' input directory needs to
#' contain the actual sharedParamNum_WASA_* files instead of the *_tpl.dat files.
#'
#' @param output_sel A named data.frame with elements 'object', 'variable', and 'digits'
#' defining the output to be generated. Needed for control file 'output_selection.txt'.
#' See the ECHSE core manual for further information. If \code{NULL} (default), the
#' variable 'out' of object 'node_su_out_1' will be given (i.e. the basin's outlet river
#' flow in m3/s).
#'
#' @param output_dbg A character vector specifying the names of objects for which
#' debug output shall be written. Needed for control file 'output_debug.txt'. If set
#' to \code{NULL} (default), no debug output will be generated.
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
#' @param error2warn Value of type \code{logical}. Shall runtime errors of the model be
#' reported as a warning instead of stopping this function with an error? Default: \code{FALSE}.
#'
#' @param nthreads Number of cores that shall be employed for the ECHSE run (argument
#' 'number_of_threads' in configuration file). See ECHSE core manual for more information.
#'
#' @return Function returns nothing.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
echse_run <- function(
  dir_input = NULL,
  run_name = NULL,
  echse_app = NULL,
  choices = NULL,
  output_sel = NULL,
  output_dbg = NULL,
  sim_start = NULL,
  sim_end = NULL,
  resolution = 86400,
  warmup_start = NULL,
  warmup_len = 3,
  max_pre_runs = 20,
  storage_tolerance = 0.01,
  error2warn = FALSE,
  nthreads = 1
) {

  # create run directory
  run_pars <- paste(dir_input, "runs", run_name, "pars", sep="/")
  run_out <- paste(dir_input, "runs", run_name, "out", sep="/")
  dir.create(paste(dir_input, "runs", run_name, sep="/"), recursive = T, showWarnings = F)
  dir.create(run_pars, recursive = T, showWarnings = F)
  dir.create(run_out, recursive = T, showWarnings = F)

  # shared parameter files path (will be overwritten if choices are given)
  sharedsvc_path <- paste(dir_input, "data", "parameter", "sharedParamNum_WASA_svc.dat", sep="/")
  sharedlu_path <- paste(dir_input, "data", "parameter", "sharedParamNum_WASA_lu.dat", sep="/")
  sharedrch_path <- paste(dir_input, "data", "parameter", "sharedParamNum_WASA_rch.dat", sep="/")



  # MODEL CHOICES #
  if(!is.null(choices)) {
    sharedsvc_path <- paste(run_pars, "sharedParamNum_WASA_svc.dat", sep="/")
    sharedlu_path <- paste(run_pars, "sharedParamNum_WASA_lu.dat", sep="/")
    sharedrch_path <- paste(run_pars, "sharedParamNum_WASA_rch.dat", sep="/")
    # SVC
    sharedpar_dat <- read.table(paste(dir_input, "data/parameter/sharedParamNum_WASA_svc_tpl.dat", sep="/"), header=T, sep="\t")

    sharedpar_dat$value <- sub("ODESOLVE", choices$odesolv, sharedpar_dat$value)
    sharedpar_dat$value <- sub("INFIL", choices$infil, sharedpar_dat$value)
    sharedpar_dat$value <- sub("EVAP", choices$evap, sharedpar_dat$value)
    sharedpar_dat$value <- sub("RCS", choices$rcs, sharedpar_dat$value)
    sharedpar_dat$value <- sub("ROUGHL", choices$roughl, sharedpar_dat$value)
    sharedpar_dat$value <- sub("DISPL", choices$displ, sharedpar_dat$value)
    sharedpar_dat$value <- sub("GLOMAX", choices$glomax, sharedpar_dat$value)
    sharedpar_dat$value <- sub("PERC", choices$perc, sharedpar_dat$value)
    sharedpar_dat$value <- sub("SOILMOD", choices$soilmod, sharedpar_dat$value)

    write.table(sharedpar_dat, sharedsvc_path,
                col.names = T, row.names = F, sep="\t", quote = F)

    # LU
    sharedparlu_dat <- read.table(paste(dir_input, "data/parameter/sharedParamNum_WASA_lu_tpl.dat", sep="/"), header=T, sep="\t")

    sharedparlu_dat$value <- sub("RUNCONC", choices$runconc, sharedparlu_dat$value)
    sharedparlu_dat$value <- sub("GROUNDWATER", choices$groundwater, sharedparlu_dat$value)

    write.table(sharedparlu_dat, sharedlu_path,
                col.names = T, row.names = F, sep="\t", quote = F)

    # RCH
    sharedparrch_dat <- read.table(paste(dir_input, "data/parameter/sharedParamNum_WASA_rch_tpl.dat", sep="/"), header=T, sep="\t")

    sharedparrch_dat$value <- sub("ROUTING", choices$routing, sharedparrch_dat$value)
    sharedparrch_dat$value <- sub("TRANSLOSS", choices$transloss, sharedparrch_dat$value)

    write.table(sharedparrch_dat, sharedrch_path,
                col.names = T, row.names = F, sep="\t", quote = F)
  } # are there choices?



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
  output_dbg <- ifelse(is.null(output_dbg), "none", output_dbg)
  content <- c("# Objects for which debug info is to be printed", "object", output_dbg)
  write(content, file = paste(run_pars, "output_debug.txt", sep="/"), sep="\n")

  # output selection file
  if(is.null(output_sel)) output_sel <- data.frame(object = "node_su_out_1", variable = "out", digits = 12)
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
      warning(paste0("ECHSE returned a runtime error during simulation, see log file: ", paste(run_out, "run.err.html", sep="/"), "."))
    } else {
      stop(paste("ECHSE returned a runtime error during simulation, see log file:", paste(run_out, "run.err.html", sep="/")))
    }
  }


} # EOF
