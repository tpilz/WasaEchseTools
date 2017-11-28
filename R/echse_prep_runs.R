#' Preparation of input for the WASA engine of the ECHSE simulation environment
#'
#' Makes use of spatial model input prepared with the \code{\link[lumpR]{lumpR}} package,
#' pre-processed meteorological data, and further input to this function to generate a
#' directory for the model ready to run simulations.
#'
#' @param sp_input_dir Character string of the directory containing the spatial input.
#' E.g. the output of \code{\link[lumpR]{db_echse_input}}.
#'
#' @param wgt_file Character string specifying the name of the file containing the
#' weighting of input locations of meteorological data. Should be the output of function
#' \code{\link[geostat]{externalInputLocationsTable}} (argument 'file_result').
#'
#' @param meteo_ext_datafiles Data.frame containing information of ECHSE's external
#' datafile input. Requires the variables 'variable', 'sums', 'past', and 'file'
#' (see the ECHSE manual for more information).
#'
#' @param prep_tpl Logical. Shall template parameter files be created? If \code{TRUE},
#'  files sharedParamNum_WASA_svc_tpl.dat, sharedParamNum_WASA_lu_tpl.dat, and
#'  sharedParamNum_WASA_rch_tpl.dat will be created (and the old files be deleted
#'  accordingly) where all choice_* parameter values will be replaced by placeholders
#'  to be adjusted, for instance, within a multi-hypothesis study employing different
#'  model structures (i.e. using different choice_* parameter realisations). Default: \code{FALSE}.
#'
#' @param echse_sim_dir Character string of the location the output shall be written
#' to (directory will be created if it does not exist).
#'
#' @return Function returns nothing.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
echse_prep_runs <- function(
  sp_input_dir = NULL,
  wgt_file = NULL,
  meteo_ext_datafiles = NULL,
  prep_tpl = FALSE,
  echse_sim_dir = getwd()) {

  # create echse simulation directory for current project
  dir.create(echse_sim_dir, recursive = T, showWarnings = F)

  # copy created echse input data
  file.copy(paste0(sp_input_dir, "/."), echse_sim_dir, recursive = T, overwrite=T)

  # dynamic variables #
  # get objDecl (created by lumpR::db_echse_input)
  objDecl_file <- paste(echse_sim_dir, "data/catchment/objDecl.dat", sep="/")

  # external location table from vegetation parameter time series (created by lumpR::db_echse_input)
  file_wgt_veg <- paste(echse_sim_dir, "data/vegPar_time_series/input_ext_locs.dat", sep="/")

  # Output file from interolation: external locations and weights table for meteo input
  file_ext_wgt <- paste(echse_sim_dir, "data/forcing/inputs_ext_locations.dat", sep="/")

  # external data locations file compiled from information given to meteo_ts_echse
  file_ext_datafiles <- paste(echse_sim_dir, "data/forcing/inputs_ext_datafiles.dat", sep="/")



  # METEO DATA #
  # prepare forcing sub-directory in echse dir
  dir.create(paste(echse_sim_dir, "data/forcing/", sep="/"))

  # Read obj declaration
  objDecl_dat <- read.table(objDecl_file, header=T)

  # read weights data
  dat_wgt <- read.table(wgt_file, header = T, sep = "\t")

  # get SVCs for each subbasin (name scheme: svc_{id_subbasin}_{id_lu}_{id_tc}_{id_svc})
  svc <- objDecl_dat$object[grep("WASA_svc", objDecl_dat$objectGroup)]

  # put SVCs into weights table replacing corresponding subbasin
  sub <- unlist(strsplit(as.character(svc), "_"))
  sub <- sub[seq(2, length(sub), 5)]
  wgt_out <- NULL
  for(s in unique(sub)) {
    r_decl <- which(sub == s)
    r_wgt <- which(dat_wgt$object == s)

    wgt_out <- rbind(wgt_out, merge(dat_wgt[r_wgt,], svc[r_decl]))
  }

  wgt_out <- wgt_out[,c("y", "variable", "location", "weight")]
  names(wgt_out)[1] <- "object"

  # read weights file of vegetation parameter time series
  wgt_veg <- read.table(file_wgt_veg, header=T)

  # combine wgt data
  wgt_out <- rbind(wgt_out, wgt_veg)

  # Create result
  write.table(x=wgt_out, file=file_ext_wgt, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

  write.table(meteo_ext_datafiles, file_ext_datafiles, col.names = T, row.names = F, sep="\t", quote=F)



  # ADJUST PARAMFUN #
  # adjust parameter functions; if there is only one line in a look-up table ECHSE will complain
  horpar_files <- list.files(paste(echse_sim_dir, "data/parameter/parFun_horpars", sep="/"))
  pos2area_files <- list.files(paste(echse_sim_dir, "data/parameter/parFun_pos2area", sep="/"))
  uh_files <- list.files(paste(echse_sim_dir, "data/parameter/parFun_uh", sep="/"))
  for (j in horpar_files) {
    dat <- read.table(paste(echse_sim_dir, "data/parameter/parFun_horpars", j, sep="/"), header=T, sep="\t")
    if(nrow(dat) > 1) {
      next
    } else {
      dat <- rbind(dat, rep(9999, ncol(dat)))
      write.table(dat, paste(echse_sim_dir, "data/parameter/parFun_horpars", j, sep="/"),
                  col.names = T, row.names = F, quote = F, sep="\t")
    }
  }
  for (j in pos2area_files) {
    dat <- read.table(paste(echse_sim_dir, "data/parameter/parFun_pos2area", j, sep="/"), header=T, sep="\t")
    if(nrow(dat) > 1) {
      next
    } else {
      dat <- rbind(dat, rep(9999, ncol(dat)))
      write.table(dat, paste(echse_sim_dir, "data/parameter/parFun_pos2area", j, sep="/"),
                  col.names = T, row.names = F, quote = F, sep="\t")
    }
  }
  for (j in uh_files) {
    dat <- read.table(paste(echse_sim_dir, "data/parameter/parFun_uh", j, sep="/"), header=T, sep="\t")
    if(nrow(dat) > 1) {
      next
    } else {
      dat <- rbind(dat, rep(9999, ncol(dat)))
      write.table(dat, paste(echse_sim_dir, "data/parameter/parFun_uh", j, sep="/"),
                  col.names = T, row.names = F, quote = F, sep="\t")
    }
  }



  if(prep_tpl) {
    # SHARED PARAMETERS #
    # SVC
    # prepare multi runs by introducing falgs into shared parameter file
    sharedpar_dat <- read.table(paste(echse_sim_dir, "data/parameter/sharedParamNum_WASA_svc.dat", sep="/"),
                                header=T, sep="\t")
    # flags into data
    sharedpar_dat$value[grep("choice_odesolve", sharedpar_dat$parameter)] <- "ODESOLVE"
    sharedpar_dat$value[grep("choice_inf", sharedpar_dat$parameter)] <- "INFIL"
    sharedpar_dat$value[grep("choice_et", sharedpar_dat$parameter)] <- "EVAP"
    sharedpar_dat$value[grep("choice_rcs", sharedpar_dat$parameter)] <- "RCS"
    sharedpar_dat$value[grep("choice_roughLen", sharedpar_dat$parameter)] <- "ROUGHL"
    sharedpar_dat$value[grep("choice_plantDispl", sharedpar_dat$parameter)] <- "DISPL"
    sharedpar_dat$value[grep("choice_gloradmax", sharedpar_dat$parameter)] <- "GLOMAX"
    sharedpar_dat$value[grep("choice_perc", sharedpar_dat$parameter)] <- "PERC"
    sharedpar_dat$value[grep("choice_soilmod", sharedpar_dat$parameter)] <- "SOILMOD"
    # write table
    write.table(sharedpar_dat, paste(echse_sim_dir, "data/parameter/sharedParamNum_WASA_svc_tpl.dat", sep="/"),
                col.names = T, row.names = F, sep="\t", quote=F)
    # remove old file
    invisible(file.remove(paste(echse_sim_dir, "data/parameter/sharedParamNum_WASA_svc.dat", sep="/")))

    # LU
    # prepare multi runs by introducing falgs into shared parameter file
    sharedparlu_dat <- read.table(paste(echse_sim_dir, "data/parameter/sharedParamNum_WASA_lu.dat", sep="/"),
                                  header=T, sep="\t")
    # flags into data
    sharedparlu_dat$value[grep("choice_runconc", sharedparlu_dat$parameter)] <- "RUNCONC"
    sharedparlu_dat$value[grep("choice_gw", sharedparlu_dat$parameter)] <- "GROUNDWATER"
    # write table
    write.table(sharedparlu_dat, paste(echse_sim_dir, "data/parameter/sharedParamNum_WASA_lu_tpl.dat", sep="/"),
                col.names = T, row.names = F, sep="\t", quote=F)
    # remove old file
    invisible(file.remove(paste(echse_sim_dir, "data/parameter/sharedParamNum_WASA_lu.dat", sep="/")))

    # RCH
    # prepare multi runs by introducing falgs into shared parameter file
    sharedparrch_dat <- read.table(paste(echse_sim_dir, "data/parameter/sharedParamNum_WASA_rch.dat", sep="/"),
                                   header=T, sep="\t")
    # flags into data
    sharedparrch_dat$value[grep("choice_route", sharedparrch_dat$parameter)] <- "ROUTING"
    sharedparrch_dat$value[grep("choice_transloss", sharedparrch_dat$parameter)] <- "TRANSLOSS"
    # write table
    write.table(sharedparrch_dat, paste(echse_sim_dir, "data/parameter/sharedParamNum_WASA_rch_tpl.dat", sep="/"),
                col.names = T, row.names = F, sep="\t", quote=F)
    # remove old file
    invisible(file.remove(paste(echse_sim_dir, "data/parameter/sharedParamNum_WASA_rch.dat", sep="/")))
  }

} # EOF
