#' Executes simulations with the WASA-SED model
#'
#' Function executes simulation with the WASA-SED model including warm-up runs to
#' start from stable model states as initial conditions.
#'
#' @param dir_run Character string giving the WASA-SED simulation directory. Needs
#' to contain sub-directories 'input' and 'output'.
#'
#' @param wasa_app Character string giving the system command of the WASA-SED application.
#'
#' @param warmup_start An object of class 'date' giving the start date of the warm-up period.
#' If \code{NULL} (default), the value in do.dat (lines 4 and 6) is used.
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
#' @param keep_warmup_states Logical value. Shall state storage files after finishing the warm-up
#' runs be stored for upcomming runs? If so, the *.stat files will be copied into
#' \code{dir_run}/input. WARNING: Existing *.stat files will be overwritten!
#' Default: \code{FALSE}.
#'
#' @param log_meta Character containg the name (and path) of a file into which information
#' about the function call will be written. See below for more information. Default:
#' \code{NULL}, i.e. no logile will be created. If the file already exists, the file
#' name to be created will be extended by a call to \code{\link{tempfile}}.
#'
#' @param keep_log Value of type \code{logical}. Shall a log file of the model run be written (to \code{dir_run})?
#' Default: \code{FALSE}.
#'
#' @param error2warn Value of type \code{logical}. Shall runtime errors of the model be
#' reported as a warning instead of stopping this function with an error? If so, the
#' model run's log file will be saved. Default: \code{FALSE}.
#'
#' @note To avoid warm-up runs, set \code{max_pre_runs} or \code{warmup_len} to zero.
#'
#' @return Function returns nothing.
#'
#' If argument \code{log_meta} was given, the logfile contains the following information:
#' \code{run_dir}: realisation of argument \code{dir_run}, \code{time_simrun}:
#' runtime of the simulation run (in secs.), \code{time_warmup}: total runtime for all
#' warm-up runs (in secs.), \code{warmup_iterations}: number of warm-up iterations,
#' \code{warmup_storchange}: relative storage change after the last warm-up iteration.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
wasa_run <- function(
  dir_run = NULL,
  wasa_app = NULL,
  warmup_start = NULL,
  warmup_len = 3,
  max_pre_runs = 20,
  storage_tolerance = 0.01,
  keep_warmup_states = FALSE,
  log_meta = NULL,
  keep_log = FALSE,
  error2warn = FALSE
) {
  if(!is.null(log_meta)) {
    f_log <- TRUE
    logfile = log_meta
    logdir <- sub(basename(logfile), "", logfile)
    if(file.exists(logfile)) logfile <- tempfile(sub(".[a-z]+$", "", basename(logfile)), tmpdir = logdir, fileext = sub("^[a-zA-Z0-9_-]+", "", basename(logfile)))
    if(logdir != "") dir.create(logdir, recursive = T, showWarnings = F)
    tryCatch(write(NULL, logfile), error = function(e) stop(paste("Could not create the given log file:", e)))
  } else {
    f_log <- FALSE
  }

  # MODIFY do.dat #
  #save original do.dat
  target_file <- paste(dir_run, "input", "do.dat", sep="/")
  file.copy(target_file,paste(target_file,".full_time",sep=""))
  # read do.dat
  do_dat <- readLines(paste(dir_run, "input", "do.dat", sep="/"))
  # save/load model states (in any case)
  do_dat[36] <- ".t. //load state of storages from files (if present) at start (optional)"
  do_dat[37] <- ".t. //save state of storages to files after simulation period (optional)"
  write.table(do_dat, file = paste(target_file,".full_time",sep=""), append = F, quote = F, row.names=F, col.names=F, sep="\t")
  # adjust according to length of warm-up period
  if(is.null(warmup_start)) {
    start_year <- as.numeric(strsplit(do_dat[4], "\t")[[1]][1])
    start_month <- as.numeric(strsplit(do_dat[6], "\t")[[1]][1])
  } else {
    start_year <- format(warmup_start, "%Y")
    start_month <- format(warmup_start, "%m")
  }
  start_date <- as.Date(paste(start_year, start_month, "01", sep="-"))
  end_date_prerun <- seq(start_date, by=paste(warmup_len-1, "month"), length=2)[2]
  do_dat[4] <- start_year
  do_dat[5] <- format(end_date_prerun, "%Y")
  do_dat[6] <- start_month
  do_dat[7] <- format(end_date_prerun, "%m")
  # re-write do.dat
  write.table(do_dat, file = target_file, append = F, quote = F, row.names=F, col.names=F, sep="\t")

  # WARM-UP RUNS #

  # conduct warm-up runs?
  time_warmup <- NA
  i_warmup <- NA
  rel_storage_change <- NA
  if(warmup_len > 0 & max_pre_runs > 0) {

    # WASA storage file
    storage_file <- paste(dir_run, "output","storage.stats",sep="/")
    # initialise water storage tracker
    if(file.exists(storage_file)) {
      storage_before <- read.table(storage_file, skip=1,header = F,row.names=1)
    } else {
      storage_before <- 0
    }
    # loop over pre-runs
    time_warmup <- system.time({
      for (i in 1:max_pre_runs) {

        # run WASA
        run_log <- system(command = paste0(wasa_app, " ", dir_run, "/input/do.dat"), intern = T)
        if(any(grepl("error", run_log, ignore.case = T))) {
          if(error2warn) {
            file_save <- paste0(tempfile("run_save_", dir_run), ".log")
            writeLines(run_log, file_save)
            warning(paste0("WASA returned a runtime error during warm-up, see log file: ", file_save, ". Continue model runs ..."))
          } else {
            writeLines(run_log, paste(dir_run, "run.log", sep="/"))
            stop(paste("WASA returned a runtime error during warm-up, see log file:", paste(dir_run, "run.log", sep="/")))
          }
        }

        # compare current water storage to storage after previous run
        storage_after <- read.table(storage_file, skip=1,header = F,row.names=1)
        rel_storage_change <- abs(sum(storage_after)-sum(storage_before))
        # avoid NaNs sum(storage_before)==0
        if (sum(storage_before)!=0) rel_storage_change <- rel_storage_change/sum(storage_before)

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
      file.copy(dir(paste(dir_run, "output", sep="/"), pattern = ".stat$|.stats$", full.names = T), paste(dir_run, "input", sep="/"), overwrite = T)
    }

  } # warmup run to be conducted?


  ### ACTUAL MODEL RUN ###
  # restore original do.dat
  file.rename(paste(target_file,".full_time",sep=""),target_file)
  # run WASA
  time_simrun <- system.time({
    run_log <- system(command = paste0(wasa_app, " ", dir_run, "/input/do.dat"), intern = T)
  })
  if(any(grepl("error", run_log, ignore.case = T))) {
    if(error2warn) {
      file_save <- paste0(tempfile("run_save_", dir_run), ".log")
      writeLines(run_log, file_save)
      warning(paste0("WASA returned a runtime error during simulation, see log file: ", file_save, "."))
    } else {
      writeLines(run_log, paste(dir_run, "run.log", sep="/"))
      stop(paste("WASA returned a runtime error during simulation, see log file:", paste(dir_run, "run.log", sep="/")))
    }
  }
  if(keep_log) writeLines(run_log, paste(dir_run, "run.log", sep="/"))

  # write logfile
  if(f_log) {
    time_end <- Sys.time()
    out_log <- data.frame(group = c(rep("meta", 5)),
                          variable = c("run_dir", "time_simrun", "time_warmup", "warmup_iterations", "warmup_storchange"),
                          value = c(dir_run, round(time_simrun["elapsed"], 1), round(time_warmup["elapsed"], 1), i_warmup, round(rel_storage_change, 3))
    )
    write.table(out_log, file=logfile, sep="\t", quote=F, row.names=F, col.names=T)
  }

} # EOF
