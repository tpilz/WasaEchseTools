#' Preparation of WASA-SED time series input by IDW interpolation of climate data to subbasin centroids
#'
#' @param file_target_locs Mandatory argument \code{file_targets} for
#' \code{\link[geostat]{externalInputLocationsTable}}. To calculate extraterrestrial
#' radiation, the file further needs to contain column 'lat' with the latitudes in
#' decimal degrees of the target locations.
#' @param file_source_locs Mandatory argument \code{file_sources} for
#' \code{\link[geostat]{externalInputLocationsTable}}.
#' @param file_result Mandatory argument for \code{\link[geostat]{externalInputLocationsTable}}.
#' File will be created within directory \code{output_dir}.
#'
#' @param file_ts Named vector of file names. Required are the four files 'precip',
#' 'temper', 'rhum', and 'glorad'. Each file is expected to hold time series of one
#' of the former variables for the various source locations (as named columns) given in
#' \code{file_source_locs}. The first column of each file, however, needs to hold a datetime
#' string in a format recognisable for function \code{\link[zoo]{read.zoo}}.
#'
#' @param output_dir Character string of the directory for \code{file_result} and the
#' WASA-SED time series input files. Will be created.
#'
#' @param resol Character string defining the resolution of the WASA-SED time series files.
#' Supported are: 'daily' (default) for daily and 'hourly' for hourly resolution. Note that
#' this affects only the time series file of precipitation and also defines the required resolution
#' of the input time series file! All other variables (temperature, humidity, radiation)
#' will be in daily resolution (and have to be given in daily resolution as well!).
#'
#' @param ... Additional arguments for \code{\link[geostat]{externalInputLocationsTable}}.
#' If not given, standard values will be used.
#'
#' @description This function performs IDW interpolation of climate data to subbasin
#' centroids (e.g. created with the \code{lumpR} package) using the
#' function \code{\link[geostat]{externalInputLocationsTable}} from the package
#' \code{\link{geostat}} which is available from the echse_tools repository: \url{https://github.com/echse/echse_tools}.
#' Extreterrestrial radiation is calculated employing function \code{\link[sirad]{extrat}}.
#'
#' @return Function returns nothing. WASA-SED time series files (and the output of
#' \code{\link[geostat]{externalInputLocationsTable}}) will be created.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
wasa_meteo <- function(
  file_target_locs = NULL,
  file_source_locs = NULL,
  file_result = NULL,
  file_ts = NULL,
  output_dir = NULL,
  resol = "daily",
  ...
) {

  # create output directory
  dir.create(output_dir, recursive = T, showWarnings = F)
  wgt_file <- paste(output_dir, file_result, sep="/")

  # calculate weights for precipitation interpolation to subbasin centroids
  externalInputLocationsTable(file_targets= file_target_locs,
                              files_sources= file_source_locs,
                              file_result= wgt_file,
                              ...)

  # read weights data
  dat_wgt <- read.table(wgt_file, header = T, sep = "\t")
  # read catchment data
  dat_cat <- read.table(file_target_locs, header = T)

  # loop over variables
  for (v in names(file_ts)) {

    # read ts data
    dat_ts <- read.zoo(file_ts[v], header=T, sep="\t", tz ="UTC")

    # remove leading "X" in colnames if present
    colnames(dat_ts) <- gsub("^X", "", colnames(dat_ts))

    # get relevant weight data
    var_wgt <- dat_wgt[grep(v, dat_wgt$variable),]

    # loop over target stations
    dat_interp <- NULL
    for (s in 1:nrow(dat_cat)) {

      # target location
      stat_tar <- dat_cat$pid[s]

      # get relevant rows in weight table
      rows_wgt <- which(var_wgt$object == stat_tar)

      # get relevant stations from ts data
      dat_t <- dat_ts[,as.character(var_wgt$location[rows_wgt]), drop=F]
      weights_s <- var_wgt$weight[rows_wgt]

      # compute weighted mean
      dat_interp_t <- apply(dat_t, 1, function(x,w=weights_s) sum(x*w))

      # combine output
      dat_interp <- cbind(dat_interp, dat_interp_t)

    } # loop over stations

    # xts object for interpolated precipitation
    dat_interp_xts <- xts(dat_interp, index(dat_ts))
    colnames(dat_interp_xts) <- dat_cat$pid

    # wasa time series input file structure
    wasa_out <- format(index(dat_interp_xts), "%d%m%Y")
    wasa_out <- cbind(wasa_out, 1:length(wasa_out))
    wasa_out <- cbind(wasa_out, round(coredata(dat_interp_xts),2))
    colnames(wasa_out) <- c("0", "0", colnames(dat_interp_xts))

    # write output
    if (v == "precip") {
      if(resol == "daily") {
        write("Daily average precipitation [mm/d] for each subasin, ordered according to Map-IDs", paste(output_dir, "rain_daily.dat", sep="/"))
        wasa_name <- "rain_daily"
      } else if(resol == "hourly") {
        write("Hourly average precipitation [mm/h] for each subasin, ordered according to Map-IDs", paste(output_dir, "rain_hourly.dat", sep="/"))
        wasa_name <- "rain_hourly"
      }
      dat_prec_interp_xts <- dat_interp_xts # needed later for rainy season
    }
    if (v == "temper") {
      write("Daily average temperature (in degree Celcius) for each subasin, ordered according to Map-IDs", paste(output_dir, "temperature.dat", sep="/"))
      wasa_name <- "temperature"
    }
    if (v == "rhum"){
      write("Daily average humidity [in %] for each subasin, ordered according to Map-IDs", paste(output_dir, "humidity.dat", sep="/"))
      wasa_name <- "humidity"
    }
    if (v == "glorad") {
      write("Daily average shortwave radiation [in W/m2] for each subasin, ordered according to Map-IDs", paste(output_dir, "radiation.dat", sep="/"))
      wasa_name <- "radiation"
    }

    write("Date  No. of days, Subasin-ID, Subasin-ID,...", paste0(output_dir, wasa_name,".dat"), append=T)
    suppressWarnings(write.table(wasa_out, paste0(output_dir, wasa_name,".dat"), col.names=T, row.names=F, append=T, quote=F, sep="\t"))

  } # loop over variables



  # Extraterrestrial radiation

  # create dummy year date vector of a non-leap year
  dummy_year <- seq.Date(as.Date("2014-01-01", format="%Y-%m-%d"), as.Date("2014-12-31", format="%Y-%m-%d"), by="day")

  # get mean latitude in dec deg.
  mean_lat <- (max(dat_cat$lat) + min(dat_cat$lat)) / 2

  # calculate extraterrestrial radiation and convert unit [MJm-2] -> [Wm-2]
  extraterr_d <- extrat(i=1:365,lat=radians(mean_lat))$ExtraTerrestrialSolarRadiationDaily * 1e6/86400

  # monthly values
  extraterr_m <- aggregate(extraterr_d, by=list(month=format(dummy_year, "%m")), mean)

  # write WASA file
  write("Extra-terrestrial shortwave radiation as monthly mean daily value in [W/m2]", paste(output_dir, "extraterrestrial_radiation.dat", sep="/"))
  write.table(extraterr_m$x, paste(output_dir, "extraterrestrial_radiation.dat", sep="/"), sep="\n", col.names=F, row.names=F, append=T)

} # EOF
