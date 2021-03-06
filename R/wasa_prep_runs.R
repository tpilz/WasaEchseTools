#' Preparation of WASA-SED model input
#'
#' Makes use of spatial model input prepared with the \code{\link[lumpR]{lumpR}} package,
#' pre-processed meteorological data, and further input to this function to generate a
#' directory for the model ready to run simulations.
#'
#' @param sp_input_dir Character string of the directory containing the spatial WASA-SED
#' input. E.g. the output of \code{\link[lumpR]{db_wasa_input}} (argument 'dest_dir').
#'
#' @param meteo_dir Character string of the directory containing the time series input
#' files for the WASA-SED model. Requires the files temperature.dat, radiation.dat,
#' humidity.dat, and rain_daily.dat or rain_hourly.dat.
#'
#' @param radex_file Character string of the file (including path) of preprared
#' extraterrestrial radiation input (named 'extraterrestrial_radiation.dat').
#'
#' @param wasa_sim_dir Character string of the location the output shall be written
#' to (directory will be created if it does not exist).
#'
#' @param proj_name Character string of the current project name, i.e. the directory
#' within \code{wasa_sim_dir} which will contain the directories 'input' (the complete
#' WASA-SED input) and 'output' (WASA-SED simulation output will be written into this
#' directory).
#'
#' @param sim_start Object of class 'date' giving the start date of the simulation
#' (will be written into WASA-SED input file 'do.dat').
#'
#' @param sim_end Object of class 'date' giving the end date of the simulation
#' (will be written into WASA-SED input file 'do.dat').
#'
#' @param timestep Integer giving the temporal resolution of the model run in hours.
#' Must be 1 (hourly) or 24 (daily resolution). Will be written into WASA-SED input
#' file 'do.dat'. Default: 24.
#'
#' @param outfiles Character vector specifying the contents of input file 'outfiles.dat'
#' controlling the output to be generated by the model. See the WASA-SED socumentation
#' for supported values. If \code{NULL} (default), the file will not be generated and
#' default output will be produced.
#'
#' @return Function returns nothing.
#'
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
wasa_prep_runs <- function(
  sp_input_dir = NULL,
  meteo_dir = NULL,
  radex_file = NULL,
  wasa_sim_dir = getwd(),
  proj_name = "wasa_run",
  sim_start = NULL,
  sim_end = NULL,
  timestep = 24,
  outfiles = NULL) {

  # adjust model input data #
  # create WASA simulation directory for current project
  dir.create(paste(wasa_sim_dir, proj_name, sep="/"), recursive = T, showWarnings = F)
  dir.create(paste(wasa_sim_dir, proj_name, "input", sep="/"), recursive = T, showWarnings = F)
  dir.create(paste(wasa_sim_dir, proj_name, "output", sep="/"), recursive = T, showWarnings = F)

  # copy created WASA input data
  file.copy(paste0(sp_input_dir, "/."), paste(wasa_sim_dir, proj_name, "input", sep="/"), recursive = T, overwrite=T)

  # do.dat
  do_dat <- readLines(paste(wasa_sim_dir, proj_name, "input", "do.dat", sep="/"))
  do_dat[2] <- "../input/"
  do_dat[3] <- "../output/"
  do_dat[4] <- paste0(format(sim_start, "%Y"), "\t//tstart (start year of simulation)")
  do_dat[5] <- paste0(format(sim_end, "%Y"), "\t//tstop (end year of simulation)")
  do_dat[6] <- paste0(format(sim_start, "%m"), "\t//mstart (start month of simulation)")
  do_dat[7] <- paste0(format(sim_end, "%m"), "\t//mstop (end month of simulation)")
  do_dat[30] <- paste0(timestep, "\t//dt: time step in [hours]")

  writeLines(do_dat, paste(wasa_sim_dir, proj_name, "input", "do.dat", sep="/"))

  # time series (meteo) data
  dir.create(paste(wasa_sim_dir, proj_name, "input", "Time_series", sep="/"), recursive = T, showWarnings = F)
  file.copy(paste0(meteo_dir, "/."), paste(wasa_sim_dir, proj_name, "input", "Time_series", sep="/"), recursive = T)
  file.copy(radex_file, paste(wasa_sim_dir, proj_name, "input", "Time_series", sep="/"))

  # ksat calibration file
  dir.create(paste(wasa_sim_dir, proj_name, "input", "Others", sep="/"), recursive = T, showWarnings = F)
  write.table(data.frame(soil_id=-1, ksat_calib_factor=1), paste(wasa_sim_dir, proj_name, "input", "Others", "calibration.dat", sep="/"),
              quote=F, row.names=F, sep="\t")

  # outfiles.dat
  content <- c("This files describe which output files are generated", outfiles)
  write(content, file = paste(wasa_sim_dir, proj_name, "input", "outfiles.dat", sep="/"), sep="\n")

  # look for existing storage files to be used
  files_stat <- dir(paste(wasa_sim_dir, proj_name, "input", sep="/"), pattern = ".stat|.stats", full.names = T)
  if(length(files_stat) > 0) file.copy(files_stat, paste(wasa_sim_dir, proj_name, "output/", sep="/"), overwrite = T)

} # EOF
