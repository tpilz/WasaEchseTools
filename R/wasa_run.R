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
#' @param warmup_len Integer giving the length of the warm-up period in months. Default: 3.
#'
#' @param max_pre_runs Integer specifying the maximum number of warm-up iterations to be
#' applied. If the relative storage change is still larger than \code{storage_tolerance}
#' after \code{max_pre_runs} iterations, the warm-up will be aborted and the model be run
#' anyway. A warning will be issued.
#'
#' @param storage_tolerance Numeric value giving the relative change of the model's water
#' storages between two connsecutive warm-up runs below which the warm-up will be
#' concluded and the actual model simulation be started.
#'
#' @return Function returns nothing.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
wasa_run <- function(
  dir_run = NULL,
  wasa_app = NULL,
  warmup_len = 3,
  max_pre_runs = 20,
  storage_tolerance = 0.01
) {

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
  start_year <- as.numeric(strsplit(do_dat[4], "\t")[[1]][1])
  start_month <- as.numeric(strsplit(do_dat[6], "\t")[[1]][1])
  start_date <- as.Date(paste(start_year, start_month, "01", sep="-"))
  end_date_prerun <- seq(start_date, by=paste(warmup_len-1, "month"), length=2)[2]
  do_dat[5] <- format(end_date_prerun, "%Y")
  do_dat[7] <- format(end_date_prerun, "%m")
  # re-write do.dat
  write.table(do_dat, file = target_file, append = F, quote = F, row.names=F, col.names=F, sep="\t")

  # WARM-UP RUNS #
  # WASA storage file
  storage_file <- paste(dir_run, "output","storage.stats",sep="/")
  # initialise water storage tracker
  storage_before <- 0
  # loop over pre-runs
  for (i in 1:max_pre_runs) {

    # run WASA
    run_log <- system(command = paste0(wasa_app, " ", dir_run, "/input/do.dat"), intern = T)
    if(any(grepl("error", run_log, ignore.case = T))) {
      writeLines(run_log, paste(dir_run, "run.log", sep="/"))
      stop(paste("WASA returned a runtime error during warm-up, see log file:", paste(dir_run, "run.log", sep="/")))
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
  if(i == max_pre_runs)
    warning(paste("Relative storage change after 'max_pre_runs' iterations was still above the tolerance threshold: ", round(rel_storage_change, 3)))

  ### ACTUAL MODEL RUN ###
  # restore original do.dat
  file.rename(paste(target_file,".full_time",sep=""),target_file)
  # run WASA
  run_log <- system(command = paste0(wasa_app, " ", dir_run, "/input/do.dat"), intern = T)
  if(any(grepl("error", run_log, ignore.case = T))) {
    writeLines(run_log, paste(dir_run, "run.log", sep="/"))
    stop(paste("WASA returned a runtime error during simulation, see log file:", paste(dir_run, "run.log", sep="/")))
  }

} # EOF
