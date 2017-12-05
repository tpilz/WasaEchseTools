#' Modify input parameters of the WASA-SED model
#'
#' Function replaces specified parameter values within prepared WASA-SED model input
#' files by given substitutes.
#'
#' @param pars A named vector of type numeric with the values of selected parameters.
#'
#' @param dir_run Character string of the path to the WASA-SED input directory.
#'
#' @details IMPORTANT: Note that some parameters are factors or summands and the
#' parameters are updated by multiplication or summation rather than mere replacement.
#' Read notes below!
#'
#' Parameter names below can be modified by adding prefix 'log_*'. In this case,
#' the value will be replaced by exp(value) internally for better handling of highly
#' uncertain parameters. For instance, 'kf_bedrock_f' might vary between 0.01 and 100
#' which can cause some numerical issues to the calibration algorithm (i.e. underrepresenting
#' small numbers < 1). In such a case it is better to specify, e.g. log(100) or log(0.01)
#' (i.e. 4.6 or -4.6, respectively) as limits for the calibration algorithm. Note that
#' only the natural logarithm (R function log()) is supported herein!
#'
#'
#' Currently supported are the following parameters (use these as names for the input
#' vector \code{pars}):
#'
#' File Hillslope/vegetation.dat (equally applied to all vegetation types and node points):
#'
#'  - stomr_f: Factor, somatal resistance\cr
#'  - rootd_f: Factor, rootdepth\cr
#'  - albedo_f: Factor, albedo
#'
#' File Hillslope/soil.dat (equally applied to all types):
#'
#'  - n_f: Factor, porosity\cr
#'  - cfr: Summand, added to coarse fragments\cr
#'
#' File River/response.dat:
#'
#'  - uhg_f: Factor, applied to both lag and response times
#'
#' File Hillslope/soter.dat:
#'
#'  - gw_delay_f: Factor, groundwater delay\cr
#'  - soildepth_f: Factor, soil depth\cr
#'  - kf_bedrock_f: Factor, kf of bedrock\cr
#'  - riverdepth_f: Factor, depth of riverbed
#'
#' File do.dat:
#'
#'  - kfcorr: Replacement of kfcorr value\cr
#'  - kfcorr0: Replacement of kfcorr0 value\cr
#'  - kfcorr_a: Replacement of kfcorr_a value\cr
#'  - intcf: Replacement of intcf value\cr
#'  - dosediment: Replacement of dosediment flag
#'
#' File calibration.dat:
#'
#'  - ksat_factor: Replacement of ksat factor (factor internally used in WASA-SED)
#'
#'
#' @return Function returns nothing.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
wasa_modify_pars <- function(
  pars = NULL,
  dir_run = NULL
) {

  # replace log values by normal values for use in WASA (log is just a matter of sampling)
  if(any(grepl(names(pars), pattern = "^log_"))) {
    log_trans = grepl(names(pars), pattern = "^log_") #find log-transformed parameters
    pars[log_trans]=exp(pars[log_trans]) #transform back to non-log scale
    names(pars) = sub(names(pars), pattern = "^log_", rep="") #remove "log_" from name
  }

  pars_t <- pars
  # parameters in vegetation.dat
  if (any(names(pars_t) %in% c("stomr_f", "rootd_f", "albedo_f"))) {
    # file that hold the parameters to be changed
    target_file <- paste(dir_run,"Hillslope/vegetation.dat",sep="/")
    # read data
    file_content <- read.table(target_file, skip=2, header = FALSE, sep = "\t", dec = ".", fill = TRUE)
    if (all(!is.finite(file_content[,ncol(file_content)]))) file_content[,ncol(file_content)]=NULL  #discard last column if empty
    # multiply somatal resistance in vegetation.dat by factor
    nn= which(names(pars_t)=="stomr_f")
    if (length(nn)>0) {
      file_content[,2]=file_content[,2]*pars_t[nn]
      pars_t <- pars_t[-nn]
    }
    # multiply rootdepth in vegetation.dat by factor (all 4 values)
    nn= which(names(pars_t)=="rootd_f")
    if (length(nn)>0) {
      file_content[,c(9:12)]=file_content[,c(9:12)]*pars_t[nn]
      pars_t <- pars_t[-nn]
    }
    # multiply albedo in vegetation.dat by factor (all 4 values)
    nn= which(names(pars_t)=="albedo_f")
    if (length(nn)>0) {
      file_content[,c(17:20)]=pmin(1,file_content[,c(17:20)]*pars_t[nn])
      pars_t <- pars_t[-nn]
    }
    #re-write file
    content=paste("Specification of vegetation parameters\nVeg-ID	Stomata_Resistance[s/m]	minsuction[hPa]	maxsuction[hPa]	height1[m]	height2[m]	height3[m]	height4[m]	rootdepth1[m]	rootdepth2[m]	rootdepth3[m]	rootdepth4[m]	LAI1[-]	LAI2[-]	LAI3[-]	LAI4[-]	albedo1[-]	albedo2[-]	albedo3[-]	albedo4[-]",sep = "")
    write(content, file = target_file)     #write header
    write.table(round(file_content, 2), file = target_file, append = TRUE, quote = F,row.names=F,col.names=F,sep="\t", na = "")
  }
  # parameters in soil.dat
  if (any(names(pars_t) %in% c("n_f", "cfr"))) {
    # file that hold the parameters to be changed
    target_file=paste(dir_run,"Hillslope/soil.dat",sep="/")
    # read data
    file_content = read.table(target_file, skip=2,header = FALSE, sep = "\t", dec = ".", fill = TRUE)
    if (all(!is.finite(file_content[,ncol(file_content)]))) file_content[,ncol(file_content)]=NULL  #discard last column if empty

    # multiply porosity in soil.dat by factor (equally for every horizon and soil)
    nn= which(names(pars_t)=="n_f")
    if (length(nn)>0) {
      # iterate over soils (lines in data)
      for (s in 1:nrow(file_content)) {
        # get number of horizons
        nhor = file_content[s,2]
        # update parameter
        file_content[s,12+(1:nhor-1)*13] = pmin(1,file_content[s,12+(1:nhor-1)*13]*pars_t[nn])
      }
      pars_t <- pars_t[-nn]
    }
    # add 'cfr' to coarse fragments in soil.dat (additive!, equally for every horizon and soil)
    nn= which(names(pars_t)=="cfr")
    if (length(nn)>0) {
      # iterate over soils (lines in data)
      for (s in 1:nrow(file_content)) {
        # get number of horizons
        nhor = file_content[s,2]
        # update parameter
        file_content[s,14+(1:nhor-1)*13] = pmin(1,pmax(0,file_content[s,14+(1:nhor-1)*13]+pars_t[nn]))
      }
      pars_t <- pars_t[-nn]
    }
    #re-write file
    content=paste("Specification of soil parameters\nSoil-ID[-]	number(horizons)[-]	(n_res[Vol-]	n_PWP[-]	n_FK2.6[-]	n_FK1.8[-]	n_nFK[-]	n_saturated[-]	n_thickness[mm]	n_ks[mm/d]	n_suction[mm]	n_pore-size-index[-]	n_bubblepressure[cm]	n_coarse_frag[-]*n	n_shrinks[0/1])	bedrock[0/1]	alluvial[0/1]",sep = "")
    write(content, file = target_file)     #write header
    write.table(round(file_content,3), file = target_file, append = TRUE, quote = F,row.names=F,col.names=F,sep="\t", na = "")
  }
  # parameters in response.dat
  if (any(names(pars_t) %in% c("uhg_f"))){
    # file that hold the parameters to be changed
    target_file=paste(dir_run,"River/response.dat",sep="/")
    # read data
    file_content = read.table(target_file, skip=2,header = FALSE, sep = "\t", dec = ".", fill = TRUE)
    if (all(!is.finite(file_content[,ncol(file_content)]))) file_content[,ncol(file_content)]=NULL  #discard last column if empty
    # multiply lag and response times in response.dat by factor
    nn= which(names(pars_t)=="uhg_f")
    if (length(nn)>0) {
      file_content[,2:3] <- file_content[,2:3]*pars_t[nn]
      pars_t <- pars_t[-nn]
    }
    #re-write file
    content=paste("Specification of routing parameter\nSubbasin-ID	lag time [d]	retention [d]",sep = "")
    write(content, file = target_file)     #write header
    write.table(round(file_content,2), file = target_file, append = TRUE, quote = F,row.names=F,col.names=F,sep="\t", na = "")
  }
  # parameters in soter.dat
  if (any(names(pars_t) %in% c("gw_delay_f","soildepth_f","kf_bedrock_f","riverdepth_f"))) {
    # file that hold the parameters to be changed
    target_file=paste(dir_run,"Hillslope/soter.dat", sep="/")
    # read data
    file_content = read.table(target_file, skip=2,header = FALSE, sep = "\t", dec = ".", fill = TRUE)
    if (all(!is.finite(file_content[,ncol(file_content)]))) file_content[,ncol(file_content)]=NULL  #discard last column if empty
    # consider shorter lines (LUs with less TCs) and bring fields to consistent position in matrix
    max_n_tcs=max(file_content[,2])
    shorter_lines=which(file_content[,2] != max_n_tcs)
    n_fields=ncol(file_content)-max_n_tcs -2 #number of fields after specification of TC-IDs
    if(length(shorter_lines) > 0) {
      for (ll in shorter_lines) {
        n_tcs = file_content[ll,2]
        file_content[ll, 2+max_n_tcs+(1:n_fields)] = file_content[ll, 2+n_tcs    +(1:n_fields)]
      }
    }
    nn= which(names(pars_t)=="gw_delay_f")
    if (length(nn)>0) {
      file_content[,ncol(file_content)]=file_content[,ncol(file_content)]*pars_t[nn]
      pars_t <- pars_t[-nn]
    }
    nn= which(names(pars_t)=="soildepth_f")
    if (length(nn)>0) {
      not_minus1=file_content[,ncol(file_content)-4]!=-1   #only modify entries not having flag=-1
      file_content[not_minus1,ncol(file_content)-4]=file_content[not_minus1,ncol(file_content)-4]*pars_t[nn]
      not_minus1=file_content[,ncol(file_content)-5]!=-1
      file_content[not_minus1,ncol(file_content)-5]=file_content[not_minus1,ncol(file_content)-5]*pars_t[nn]
      pars_t <- pars_t[-nn]
    }
    nn= which(names(pars_t)=="kf_bedrock_f")
    if (length(nn)>0) {
      file_content[,ncol(file_content)-7]=file_content[,ncol(file_content)-7]*pars_t[nn]
      pars_t <- pars_t[-nn]
    }
    nn= which(names(pars_t)=="riverdepth_f")
    if (length(nn)>0) {
      file_content[,ncol(file_content)-3]=file_content[,ncol(file_content)-3]*pars_t[nn]
      pars_t <- pars_t[-nn]
    }
    #consider shorter lines (LUs with less TCs) and bring fields to "spars_te" representation (unequal number of fields)
    if(length(shorter_lines) > 0) {
      for (ll in shorter_lines) {
        file_content[ll,]
        n_tcs = file_content[ll,2]
        file_content[ll, 2+n_tcs    +(1:n_fields)] =
          file_content[ll, 2+max_n_tcs+(1:n_fields)]
        file_content[ll, ncol(file_content)+1-1:(max_n_tcs-n_tcs)] = NA #mask obsolete fields with NA
      }
    }
    #re-write file
    content=paste("Specification of landscape units\nLU-ID[id]  No._of_TC[-]	TC1[id]	TC2[id]	TC3[id]	kfsu[mm/d]	length[m]	meandep[mm]	maxdep[mm]	riverbed[mm]	gwflag[0/1]	gw_dist[mm]	frgw_delay[day]",sep = "")
    write(content, file = target_file)     #write header
    write.table(round(file_content, 1), file = target_file, append = TRUE, quote = F,row.names=F,col.names=F,sep="\t", na = "")
  }
  # params in do.dat
  if (any(names(pars_t) %in% c("kfcorr","kfcorr0","kfcorr_a", "intcf", "dosediment"))) {
    # adjust kfcorr in do.dat
    target_file=paste(dir_run,"do.dat",sep="/") #file that hold the parameters to be changed
    # read data
    file_content = read.table(target_file, skip=0,header = FALSE, sep = "$",stringsAsFactors=FALSE)
    nn= which(names(pars_t)=="kfcorr")
    if (length(nn)>0) {
      file_content[24,1]=paste(round(pars_t[nn],2),"  //kfcorr:  hydraulic conductivity factor (for daily model version) (kfcorr)",sep="")
      pars_t <- pars_t[-nn]
    }
    if(any(names(pars_t)=="kfcorr0") || any(names(pars_t)=="kfcorr_a")) {
      nn= which(names(pars_t)=="kfcorr0")
      if (length(nn)>0) {
        file_content[24,1]=paste(round(pars_t[nn],2)," ",sep="")
        pars_t <- pars_t[-nn]
      } else { # default value
        file_content[24,1]="1 "
      }
      nn= which(names(pars_t)=="kfcorr_a")
      if (length(nn)>0) {
        file_content[24,1]=paste(file_content[24,1],round(pars_t[nn],2)," 0	//kfcorr:  hydraulic conductivity factor (for daily model version) (kfcorr0) [optional: a <tab> b for kfcorr=kfcorr0*(a*1/daily_precip+b+1) ",sep="")
        pars_t <- pars_t[-nn]
      } else { # default value
        file_content[24,1]=paste(file_content[24,1],"0 0	//kfcorr:  hydraulic conductivity factor (for daily model version) (kfcorr0) [optional: a <tab> b for kfcorr=kfcorr0*(a*1/daily_precip+b+1) ",sep="")
      }
    }
    nn= which(names(pars_t)=="intcf")
    if (length(nn)>0) {
      file_content[25,1]=paste(round(pars_t[nn], 2)," //intcf: interception capacity per unit LAI (mm)",sep="")
      pars_t <- pars_t[-nn]
    }
    nn= which(names(pars_t)=="dosediment")
    if (length(nn)>0) {
      file_content[31,1]=paste(pars_t[nn],"  //dosediment")
      pars_t <- pars_t[-nn]
    }
    #rewrite file
    write.table(file_content, file = target_file, append = F, quote = F,row.names=F,col.names=F,sep="\t")
  }
  # params in calibration.dat
  if (any(names(pars_t) %in% c("ksat_factor"))) {
    # file that holds the parameters to be changed
    target_file=paste(dir_run,"Others/calibration.dat",sep="/")
    # read data
    file_content = read.table(target_file, skip=1,row.names=NULL, header = FALSE, sep = "\t", dec = ".")
    nn= which(names(pars_t)=="ksat_factor")
    if (length(nn)>0) {
      file_content[,2]=pars_t[nn]
      pars_t <- pars_t[-nn]
    }
    content=paste("Soil-ID","Ksat-calib-factor",sep="\t")
    write(content, file = target_file)     #write header
    write.table(round(file_content,3), file = target_file, append = TRUE, quote = F,row.names=F,col.names=F,sep="\t")
  }
  # all parameters used?
  if(length(pars_t) > 0)
    stop(paste("The following parameter(s) was/were not used (no instruction implemented):", paste(names(pars_t), collapse = ", ")))

} # EOF
