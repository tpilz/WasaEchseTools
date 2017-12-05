#' Comparison of simulated discharge of the WASA-SED model and the WASA engine of the ECHSE environment with observations
#'
#' Function creates a plot as pdf file showing the simulated vs. observed discharge time series.
#'
#' @return Function returns nothing. A plot (pdf file) will be created.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
plot_discharge <- function(
  file_obs = NULL,
  file_echse = NULL,
  file_wasa = NULL,
  var_obs = NULL,
  var_echse = NULL,
  var_wasa = NULL,
  dir_out = NULL,
  plot_name = NULL
) {

  # read data
  if(!is.null(file_obs)) {
    dat_obs <- read.table(file_obs, header=T, sep="\t") %>%
      mutate(date = as.POSIXct(.[[1]], tz ="UTC"), group = "obs") %>%
      rename(value = !!var_obs) %>%
      select(date, group, value)
  } else dat_obs <- NULL
  if(!is.null(file_echse)) {
    dat_echse <- read.table(file_echse, header=T, sep="\t") %>%
      mutate(date = as.POSIXct(.[[1]], tz ="UTC"), group = "echse") %>%
      rename(value = !!var_echse) %>%
      select(date, group, value)
  } else dat_echse <- NULL
  if(!is.null(file_wasa)) {
    dat_wasa <- read.table(file_wasa, header=T, skip=1, check.names = F) %>%
      mutate(date = as.POSIXct(paste(year, day, sep="-"), "%Y-%j", tz ="UTC"), group = "wasa") %>%
      rename(value = !!var_wasa) %>%
      select(date, group, value)
  } else dat_wasa <- NULL

  # combine and common time frame
  dat_all <- rbind(dat_obs, dat_echse, dat_wasa) %>%
    # replace -9999 by NA
    mutate(value = replace(value, value == -9999, NA)) %>%
    # calculate min and max common dates and combine to main data.frame()
    cbind(group_by(., group) %>%
    summarise(date_min = min(date), date_max = max(date)) %>%
    ungroup() %>%
    summarise(date_min = max(date_min), date_max = min(date_max))) %>%
    # filter common dates
    filter(date >= max(date_min) & date <= min(date_max))

  # plot
  gp <- ggplot(dat_all, aes(x=date, y = value, group=group, colour=group)) +
    geom_line() +
    scale_x_datetime(date_minor_breaks = "1 week", date_breaks = "1 month", date_labels = "%b") +
    theme_bw()

  if(is.null(plot_name)) plot_name <- "plot_discharge"
  ggsave(paste0(dir_out, "/", plot_name, ".pdf"), gp, height=8, width=14)

} # EOF
