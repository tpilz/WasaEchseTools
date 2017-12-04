#' Comparison of simulation outputs of the WASA-SED model and the WASA engine of the ECHSE environment
#'
#' Function creates a plot as pdf file showing the simulated time series of the most important
#' water balance components.
#'
#' @return Function returns nothing. A plot (pdf file) will be created.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
compare_wasa_echse <- function(
  resdir_echse = NULL,
  resdir_wasa = NULL,
  dir_out = NULL,
  plot_name = NULL,
  sub_pars = NULL,
  prec_file_wasa = NULL,
  resol = "hourly"
) {

  # determine conversion factor from x/s to x/timestep
  if(resol == "daily") {
    steplen <- 86400
  } else if(resol == "hourly") {
    steplen <- 3600
  }

  # get precipitation data (should be the same for echse and wasa)
  dat_prec <- left_join(read.table(prec_file_wasa, header=T, skip=2, check.names = F)[,-2] %>%
    mutate(date = as.POSIXct(sprintf(.[[1]], fmt="%08d"), "%d%m%Y", tz="UTC")) %>%
    select(-1) %>%
    melt(id.vars="date", variable.name = "object") %>%
    mutate(object = paste0("sub_", object), variable = "precip"),
    sub_pars,
    by = "object") %>%
    # in case of hourly resolution, get hour of day (not in file)
    group_by(date, variable, object) %>%
    mutate(hour = 0:(n()-1)) %>%
    ungroup() %>%
    mutate(date = as.POSIXct(paste0(date, " ", hour, ":00:00"), tz="UTC")) %>%
    group_by(date, variable) %>%
    # calculate catchment area
    left_join(.,
              summarise(., area_sum = sum(area)),
              by = c("variable", "date")) %>%
    # area-weighted sums over the catchment
    summarise(value = sum(value * area/area_sum))



  # aggregate ECHSE values
  dat_echse <- sub_pars %>%
    # read data, select relevant balance variables
    ddply("object", function(x) {
      read.table(paste0(resdir_echse, "/", x$object, ".txt"), header=T, sep="\t") %>%
        select(end_of_interval, etp, eta, eti, r_out_surf, r_out_inter, r_out_base, run_gw, v_soilwat) %>%
        mutate(area = x$area, eta=eta+eti, run_surf = r_out_surf, run_sub = r_out_inter + r_out_base, gw_rchrg = run_gw,
               date = as.POSIXct(end_of_interval, tz="UTC")-steplen) %>% # convert date to "begin of interval" (as in WASA output)
        select(-end_of_interval, -eti, -r_out_surf, -r_out_inter, -r_out_base, -run_gw)
    }) %>%
    # tidy data.frame()
    melt(id.vars = c("object", "area", "date"))  %>%
    group_by(date, variable) %>%
    left_join(.,
              summarise(., area_sum = sum(area)),
              by = c("variable", "date")) %>%
    # convert all values into mm/timestep for every catchment
    mutate(value = if_else(grepl("run_", variable), value*1000*steplen/(area*1e6), if_else(grepl("v_", variable), value*1000, value*1000*steplen))) %>%
    # calculate catchment-wide area-weighted sums (mm/timestep)
    summarise(value = sum(value * area/area_sum))


  # agregate WASA values
  wasa_files <- data.frame(variable = c("etp", "eta", "run_surf", "run_sub", "gw_rchrg", "v_soilwat"),
                           file = c("potetranspiration.out", "actetranspiration.out", "total_overlandflow.out",
                                    "subsurface_runoff.out", "deep_gw_recharge.out", "daily_theta.out"))
  dat_wasa <- left_join(wasa_files %>%
    # read data
    ddply("file", function(x) {
      read.table(paste(resdir_wasa, x$file, sep="/"), header=T, skip=1, check.names = F) %>%
        mutate(date = as.POSIXct(paste(Year, .[[2]], sep="-"),  "%Y-%j", tz="UTC")) %>%
        select(-matches("year|day|timestep")) %>%
        melt(id.vars="date", variable.name = "object") %>%
        mutate(object = paste0("sub_", object), variable = x$variable)
    }) %>%
    select(-file) %>%
    # in case of hourly resolution, get hour of day (not in file)
    group_by(date, variable, object) %>%
    mutate(hour = 0:(n()-1)) %>%
    ungroup() %>%
    mutate(date = as.POSIXct(paste0(date, " ", hour, ":00:00"), tz="UTC")),
    # merge with subbasin areas
    sub_pars,
    by = "object") %>%
    group_by(date, variable) %>%
    # calculate catchment area
    left_join(.,
              summarise(., area_sum = sum(area)),
              by = c("variable", "date")) %>%
    # all values in mm/timestep
    mutate(value = if_else(grepl("run_|gw_", variable), value*1000/(area*1e6), value)) %>%
    # area-weighted sums over the catchment
    summarise(value = sum(value * area/area_sum))


  # combine
  dat_all <- rbind(dat_wasa %>%
                        mutate(model = "wasa"),
                   dat_echse %>%
                         mutate(model = "echse"),
                   dat_prec %>%
                     filter(date %in% dat_echse$date) %>%
                     mutate(model = "both"))

  # summary statistics
  dat_text <- dat_all %>%
    group_by(model, variable) %>%
    mutate(value2 = if_else(grepl("v_", variable), c(diff(value), 0), value)) %>%
    summarise(x = min(date), y = 0.8*max(value), label = sum(value2)) %>%
    group_by(variable) %>%
    summarise(x = min(x), y = max(y), label = paste("Sums:", paste(model, round(label, 2), sep=" = ", collapse = ", ")))

  # # water balance
  # gw_stor <- read.table(paste(resdir_wasa, "gw_storage.stat", sep="/"), header=T, skip=1) %>%
  #   mutate(value = volume * area/sum(area))
  # gw_stor <- sum(gw_stor$value)
  # gw_stor_a <- read.table(paste(resdir_wasa, "gw_storage.stat_start", sep="/"), header=T, skip=1) %>%
  #   mutate(value = volume * area/sum(area))
  # gw_stor_a <- sum(gw_stor_a$value)
  # gw_change_wasa <- gw_stor-gw_stor_a
  #
  # dat_summary <- dat_all %>%
  #   filter(variable != "etp") %>%
  #   group_by(model, variable) %>%
  #   mutate(value2 = if_else(grepl("v_", variable), c(diff(value), 0), value)) %>%
  #   group_by(model) %>%
  #   summarise(value = sum(value2)) %>%
  #   mutate(value = if_else(grepl("precip", variable), -1*value, value)) %>%
  #   summarise(value = sum(value))

  # plot
  gp <- ggplot(data = dat_all, mapping = aes(x = date, y = value, group = model, colour = model)) +
    geom_line() +
    scale_x_datetime(date_minor_breaks = "1 week", date_breaks = "1 month", date_labels = "%b") +
    theme_bw() +
    facet_grid(variable ~ ., scales = "free_y") +
    geom_text(data = dat_text, mapping = aes(x = x, y = y, label = label, hjust=0), inherit.aes=FALSE)

  if(is.null(plot_name)) plot_name <- "wasa_echse_compare"
  ggsave(paste0(dir_out, "/", plot_name, ".pdf"), gp, height=10, width=13)


} # EOF
