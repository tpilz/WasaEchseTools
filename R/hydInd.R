#' Calculation of hydrological indices / streamflow signatures
#'
#' @param dat.q.ts Object of class 'xts' containing a discharge time series. Unit
#' must be (m3/s) whereas the resolution must be the same as for dat.pr.ts.
#'
#' @param dat.pr.ts Object of class 'xts' containing a rainfall time series. Unit
#' must be (m3) whereas the resolution must be the same as for dat.q.ts.
#'
#' @param na.rm Object of type \code{logical}. Shall \code{NA} values be removed from
#' \code{dat.q.ts} and \code{dat.pr.ts}? If \code{FALSE}, an error might occur for
#' certain calculations if \code{NA}s are found. Default: \code{TRUE}.
#'
#' @param ignore.zeros Object of type \code{logical}. Shall value less than or equal
#' to \code{thresh.zero} be removed from \code{dat.q.ts}? Will only be applied for
#' certain quantile-based calculations as many small / zero values may distort the
#' result.
#'
#' @param thresh.zero Values of \code{dat.q.ts} below this value will be treated as
#' zero flows (directed to \code{\link[hydrostats]{high.spells}} argument \code{ctf.threshold}).
#'
#' @param flood.thresh Numeric value giving the threshold in (V/T) for the definition of
#' a flood event (directed to \code{\link[hydrostats]{high.spells}} argument \code{threshold}).
#'
#' @details Function employs the package \code{\link[hydrostats]{hydrostats}} for
#' certain calculations.
#'
#' @return A named vector of type numeric with the following elements:
#'
#' rc: Runoff coefficient in (\%)
#'
#' m_h: Average annual maximum flow in (m3/s)
#'
#' m_l: Average annual minimum flow in (m3/s)
#'
#' f_h: Average number of high spell events per year
#'
#' f_l: Average number of low flow events per year
#'
#' r_r: Average absolute daily change during rise periods within high flow events (m3/s/day)
#'
#' r_f: Average absolute daily change during falling periods within high flow events (m3/s/day)
#'
#' fdc_s: Average slope of flow duration curve for medium range (33% to 66%) of non-zero flows;
#' the higher the value the more variable the flow regime, see Sawicz et al. (2011)
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
hydInd <- function(
  dat.q.ts = NULL,
  dat.pr.ts = NULL,
  na.rm = T,
  ignore.zeros = F,
  thresh.zero = NULL,
  flood.thresh = NULL
) {

  # argument checks
  if(!any(grepl("xts", class(dat.q.ts))))
    stop("Argument dat.q.ts must be an xts object!")
  if(!any(grepl("xts", class(dat.pr.ts))))
    stop("Argument dat.pr.ts must be an xts object!")
  if(any(index(dat.pr.ts) != index(dat.q.ts)))
    stop("Dates (index) of arguments dat.q.ts and dat.pr.ts must be exactly the same!")
  if(length(unique(diff(index(dat.q.ts)))) != 1 | length(unique(diff(index(dat.pr.ts)))) != 1)
    stop("Resolutions of dat.q.ts and/or dat.pr.ts are not unique!")

  # remove NAs (if TRUE)
  if(na.rm && any(is.na(dat.q.ts))) {
    nas <- which(is.na(dat.q.ts))
    dat.q.ts <- dat.q.ts[-nas]
    dat.pr.ts <- dat.pr.ts[-nas]
  }

  # remove zeros if desired (only for certain quantile-based calculations)
  if(ignore.zeros) {
    dat.q.ts.nozero <- dat.q.ts[!(dat.q.ts <= thresh.zero)]
  } else {
    dat.q.ts.nozero <- dat.q.ts
  }

  # get resolution in seconds
  resol <- as.numeric(difftime(index(dat.q.ts)[2], index(dat.q.ts)[1], units="secs"))

  # data.frame object required for external hydrostats-functions
  df.hydrostats <- data.frame(Q=as.numeric(dat.q.ts), Date=index(dat.q.ts))
  # calculate hydrostats-based statistics (relevant information assigned later on)
  low.hydrostats <- low.spells(df.hydrostats, plot = F)
  high.hydrostats <- high.spells(df.hydrostats, threshold=flood.thresh, ignore.zeros = ignore.zeros, ctf.threshold = thresh.zero, plot = F)

  # MAGNITUDE
  # runoff ratio over the whole given period [%]
  # calculate runoff ratio over the whole period
  rc <- sum(dat.q.ts)*resol / sum(dat.pr.ts) * 100
  # low flow: average annual minimum flow [V/T]
  m_l <- low.hydrostats$avg.min.ann
  # high flow: average annual maximum flow [V/T]
  m_h <- high.hydrostats$avg.max.ann

  # FLOW REGIME
  # average slope of flow duration curve for medium range (33% to 66%) of non-zero flows
  # the higher the value the more variable the flow regime; see Sawicz et al. (2011)
  quants <- quantile(dat.q.ts.nozero*100, probs=c(1-0.33,1-0.66))
  fdc_s <- (log(quants[1]) - log(quants[2])) / (0.66-0.33)

  # FREQUENCY
  # low flow: average number of low flow events per year [1/year]
  f_l <- low.hydrostats$low.spell.freq
  # high flow: average number of high spell events per year [1/year]
  f_h <- high.hydrostats$spell.freq

  # CONCENTRATION
  # Rate of change in flow events
  # averge absolute daily change during rise periods within high flow events [V/T/day]
  r_r <- high.hydrostats$avg.rise
  r_f <- high.hydrostats$avg.fall

  # output
  out <- c(rc, m_h, m_l, f_h, f_l, r_r, r_f, fdc_s)
  names(out) <- c("rc", "m_h", "m_l", "f_h", "f_l", "r_r", "r_f", "fdc_s")
  return(out)

} # EOF
