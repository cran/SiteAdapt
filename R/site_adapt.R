#'
#'  Site Adaptation of Solar Irradiance Modeled Series with Coincident Ground Measurements
#'
#' @param Target Dataframe object with solar radiation series to be adapted including time (with same time zone as subset_calibrating_period), the solar irradiance modeled series to be site adapted, along with their clear sky index and solar elevation (degrees)
#' @param latitude_target Site latitude of solar radiation series to be adapted (degrees, +N)
#' @param z_target Site elevation above sea level of solar radiation series to be adapted (m)
#' @param Calibration Dataframe object with solar radiation series for calibrating including time (with same time zone as subset_target_period), solar irradiance modeled and measured series, along with modeled clear sky index and solar elevation (degrees)
#' @param latitude_calibrat Site latitude of solar radiation series for calibrating (degrees, +N)
#' @param z_calibrat Site elevation above sea level of solar radiation series for calibrating (m)
#' @param timezone Time zone specification of the calibration_period and target_period datasets
#' @param GHI_threshold Upper limit of GHI series (same units that Target). For automatic calulation from observed data, set it to -99
#' @param DNI_threshold Upper limit of DNI series (same units that Target). For automatic calulation from observed data, set it to -99
#'
#' @return Dataframe object including time and site adapted solar irradiance series
#'
#' @examples
#' # A site located in the the Namib Desert of Namibia (Gobabeb, GOB) is selected
#'
#' # - latitude:   23.5614 S
#' # - Longitude:  15.0420 E
#' # - Elevation:  407.0 m asl
#'
#'  # Load calibration and modeled datasets
#'    data(calibration_2016) # Measured from BSRN
#'    data(target_2013_2016) # Provided by CAMS-RAD service
#'
#' observed_2013_2016$time = as.POSIXct(
#' paste(observed_2013_2016$Year, "-",
#' observed_2013_2016$Month, "-",
#' observed_2013_2016$Day, " ",
#' observed_2013_2016$Hour, ":",
#' observed_2013_2016$Minute, sep=""),
#' tz ="UTC")
#'
#' target_2013_2016$time = as.POSIXct(
#' paste(target_2013_2016$Year, "-",
#' target_2013_2016$Month, "-",
#' target_2013_2016$Day, " ",
#' target_2013_2016$Hour, ":",
#' target_2013_2016$Minute, sep=""),
#' tz ="UTC")
#'
#' \donttest{
#'
#'  # Apply the site adaptation procedure
#'    site_adapted_series = site_adapt(
#'    Target = target_2013_2016,
#'    latitude_target = -23.5614, # Latitude of target site
#'    z_target = 407.0, # Elevation of target site
#'    Calibration = calibration_2016,
#'    latitude_calibrat = -23.5614, # Same location of target period
#'    z_calibrat = 407.0, # Same location of target period
#'    timezone = "UTC",
#'    GHI_threshold = -99, # The threshold is calculated from observed data
#'    DNI_threshold = -99) # The threshold is calculated from observed data
#'
#'
#' # Load measured data for evaluating the site adaptation performance:
#'    data(observed_2013_2016)
#'
#' # Merge datasets
#'
#' site_adapted_series$time = as.POSIXct(
#' paste(site_adapted_series$Year, "-",
#' site_adapted_series$Month, "-",
#' site_adapted_series$Day, " ",
#' site_adapted_series$Hour, ":",
#' site_adapted_series$Minute, sep=""),
#' tz ="UTC")
#'
#' meas_model = merge(observed_2013_2016[,6:9],
#' target_2013_2016[,c(6:9,11)],
#' by = "time", all = FALSE )
#'
#' meas_model_adapt = merge(meas_model,
#' site_adapted_series[,6:10],
#' by = "time", all = FALSE )
#'
#'
#' # Display scatterplots
#' library(RColorBrewer)
#' pal <- rev(brewer.pal(11,"Spectral"))
#' pal=pal[2:11]
#'
#' library(ggplot2)
#' scatter_DNI.obs = ggplot() +
#' geom_hex(data=meas_model_adapt,aes(x=DNI.obs, y = DNI.mod),bins = 125, alpha = 1) +
#' scale_fill_gradientn(colours = pal)+ theme_light() +
#' xlab(expression(paste("Measured DNI (W / ", m^2, " )", sep=""))) +
#'  ylab(expression(paste("Modeled DNI (W / ", m^2, " )", sep=""))) +
#'  theme(legend.position = "none") +
#'  xlim(100, 1120) + ylim(100,1120) +
#' geom_abline(intercept = 0, slope = 1, color="purple", linetype="solid", size=1.5, alpha = 0.5)
#'
#' plot_DNI_adapt = ggplot() +
#' geom_hex(data=meas_model_adapt,aes(x=DNI.obs, y = DNI_adapted),bins = 125, alpha = 1) +
#'  scale_fill_gradientn(colours = pal)+ theme_light() +
#'  xlab(expression(paste("Measured DNI (W / ", m^2, " )", sep=""))) +
#'  ylab(expression(paste("Adapted DNI (W / ", m^2, " )", sep=""))) +
#'  theme(legend.position = "none") + xlim(100, 1120) + ylim(100,1120) +
#'  geom_abline(intercept = 0, slope = 1, color="purple", linetype="solid", size=1.5, alpha = 0.5)
#'
#' library(ggpubr)
#' ggarrange(scatter_DNI.obs, plot_DNI_adapt)
#'
#' scatter_GHI.obs = ggplot() +
#' geom_hex(data=meas_model_adapt,aes(x=GHI.obs, y = GHI.mod),bins = 125, alpha = 1) +
#'  scale_fill_gradientn(colours = pal)+ theme_light() +
#'  xlab(expression(paste("Measured GHI (W / ", m^2, " )", sep=""))) +
#'  ylab(expression(paste("Modeled GHI (W / ", m^2, " )", sep=""))) +
#'  theme(legend.position = "none") + xlim(100, 1180) + ylim(100,1180) +
#'  geom_abline(intercept = 0, slope = 1, color="purple", linetype="solid", size=1.5, alpha = 0.5)
#'
#' plot_GHI_adapt = ggplot() +
#' geom_hex(data=meas_model_adapt,aes(x=GHI.obs, y = GHI_adapted),bins = 125, alpha = 1) +
#'  scale_fill_gradientn(colours = pal) + theme_light() +
#'  xlab(expression(paste("Measured GHI (W / ", m^2, " )", sep=""))) +
#'  ylab(expression(paste("Adapted GHI (W / ", m^2, " )", sep=""))) +
#'  theme(legend.position = "none") + xlim(100, 1180) + ylim(100,1180) +
#'  geom_abline(intercept = 0, slope = 1, color="purple", linetype="solid", size=1.5, alpha = 0.5)
#' ggarrange(scatter_GHI.obs, plot_GHI_adapt)
#'
#' scatter_DHI.obs = ggplot() +
#' geom_hex(data=meas_model_adapt,aes(x=DHI.obs, y = DHI.mod),bins = 125, alpha = 1) +
#'  scale_fill_gradientn(colours = pal)+ theme_light() +
#'  xlab(expression(paste("Measured DHI (W / ", m^2, " )", sep=""))) +
#'  ylab(expression(paste("Modeled DHI (W / ", m^2, " )", sep=""))) +
#'  theme(legend.position = "none") + xlim(25, 700) + ylim(25, 700) +
#'  geom_abline(intercept = 0, slope = 1, color="purple", linetype="solid", size=1.5, alpha = 0.5)
#'
#' plot_DHI_adapt = ggplot() +
#' geom_hex(data=meas_model_adapt,aes(x=DHI.obs, y = DHI_adapted),bins = 125, alpha = 1) +
#'  scale_fill_gradientn(colours = pal)+ theme_light() +
#'  xlab(expression(paste("Measured DHI (W / ", m^2, " )", sep=""))) +
#'  ylab(expression(paste("Adapted DHI (W / ", m^2, " )", sep=""))) +
#'  theme(legend.position = "none") + xlim(25, 700) + ylim(25, 700) +
#'  geom_abline(intercept = 0, slope = 1, color="purple", linetype="solid", size=1.5, alpha = 0.5)
#' ggarrange(scatter_DHI.obs, plot_DHI_adapt)
#'
#'
#'
#'
#'  # Display ECDF plots
#' plot_ECDF_DNI = ggplot(data=meas_model_adapt[which(meas_model_adapt$Elev > 0),])+
#'  stat_ecdf(aes(DNI.obs), col="firebrick", lwd = 0.75) +
#'  stat_ecdf(aes(DNI.mod), col="dodgerblue", lwd = 0.75) +
#'  stat_ecdf(aes(DNI_adapted), col="purple", lwd = 0.75) +
#'  theme_light() + xlab(expression(paste("DNI (W / ", m^2, " )", sep="")))+ ylab("ECDF ( - )")+
#'  annotate("text", x = 50, y = 0.9, label = "Measured", col="firebrick1", size = 4)+
#'  annotate("text", x = 50, y = 0.8, label = "Modeled", col="dodgerblue", size = 4)+
#'  annotate("text", x = 50, y = 0.7, label = "Adapted", col="purple", size = 4)
#' plot_ECDF_DNI
#'
#' plot_ECDF_GHI = ggplot(data=meas_model_adapt[which(meas_model_adapt$Elev > 0),])+
#'  stat_ecdf(aes(GHI.obs), col="firebrick", lwd = 0.75) +
#'  stat_ecdf(aes(GHI.mod), col="dodgerblue", lwd = 0.75) +
#'  stat_ecdf(aes(GHI_adapted), col="purple", lwd = 0.75) +
#'  theme_light() + xlab(expression(paste("GHI (W / ", m^2, " )", sep="")))+ ylab("ECDF ( - )")+
#'  annotate("text", x = 50, y = 0.9, label = "Measured", col="firebrick1", size = 4)+
#'  annotate("text", x = 50, y = 0.8, label = "Modeled", col="dodgerblue", size = 4)+
#'  annotate("text", x = 50, y = 0.7, label = "Adapted", col="purple", size = 4)
#' plot_ECDF_GHI
#'
#' plot_ECDF_DHI = ggplot(data=meas_model_adapt[which(meas_model_adapt$Elev > 0),])+
#'  stat_ecdf(aes(DHI.obs), col="firebrick", lwd = 0.75) +
#'  stat_ecdf(aes(DHI.mod), col="dodgerblue", lwd = 0.75) +
#'  stat_ecdf(aes(DHI_adapted), col="purple", lwd = 0.75) +
#'  theme_light() + xlab(expression(paste("DHI (W / ", m^2, " )", sep="")))+ylab("ECDF ( - )")+
#'  annotate("text", x = 25, y = 0.9, label = "Measured", col="firebrick1", size = 4)+
#'  annotate("text", x = 25, y = 0.8, label = "Modeled", col="dodgerblue", size = 4)+
#'  annotate("text", x = 25, y = 0.7, label = "Adapted", col="purple", size = 4)
#' plot_ECDF_DHI
#'
#'
#'
#'
#' # Statistical indicators
#' library(hydroGOF)
#' pbias(meas_model_adapt$GHI.mod,meas_model_adapt$GHI.obs)
#' pbias(meas_model_adapt$GHI_adapted,meas_model_adapt$GHI.obs)
#'
#' pbias(meas_model_adapt$DNI.mod,meas_model_adapt$DNI.obs)
#' pbias(meas_model_adapt$DNI_adapted,meas_model_adapt$DNI.obs)
#'
#' pbias(meas_model_adapt$DHI.mod,meas_model_adapt$DHI.obs)
#' pbias(meas_model_adapt$DHI_adapted,meas_model_adapt$DHI.obs)
#'
#' rmse(meas_model_adapt$GHI.mod,meas_model_adapt$GHI.obs)
#' rmse(meas_model_adapt$GHI_adapted,meas_model_adapt$GHI.obs)
#'
#' rmse(meas_model_adapt$DNI.mod,meas_model_adapt$DNI.obs)
#' rmse(meas_model_adapt$DNI_adapted,meas_model_adapt$DNI.obs)
#'
#' rmse(meas_model_adapt$DHI.mod,meas_model_adapt$DHI.obs)
#' rmse(meas_model_adapt$DHI_adapted,meas_model_adapt$DHI.obs)}
#'
#'
#'
#' @export
#' @import glmulti
#' @import solaR
#' @import hyfo
#' @import hydroGOF
#' @import ggplot2
#' @import RColorBrewer
#' @import ggpubr
#' @import stats



site_adapt <- function(Target,latitude_target,z_target,
                       Calibration,latitude_calibrat,z_calibrat,
                       timezone,
                       GHI_threshold,
                       DNI_threshold){

# Values above the solar elevation threshold
  threshold = 2
  calibrating_period = Calibration[which(Calibration$Elev > threshold),]
  target_period = Target[which(Target$Elev > threshold),]

# Definition of turbidity thresholds
  threshold_aod_380 = 0.11
  threshold_aod_500 = 0.15

# Separation of clear and cloud conditions for the target period
  target_period$time = as.POSIXct(paste(target_period$Year, "-",
                                target_period$Month, "-",
                                target_period$Day, " ",
                                target_period$Hour, ":",
                                target_period$Minute, sep=""),
                                tz = timezone)
  top = TOA(latitude_target,target_period$Elev,target_period$time)
  clear_sky_target_DNI = SOLIS(top,target_period$Elev,
                                 aod_380 = rep(threshold_aod_380, length(target_period$Minute)),
                                 aod_500 = rep(threshold_aod_500, length(target_period$Minute)),
                                 w = rep(1, length(target_period$Minute)),
                                 latitude_target)
  clear_data_target = which(target_period$DNI.mod >= clear_sky_target_DNI)
  cloud_data_target = which(target_period$DNI.mod < clear_sky_target_DNI)

# Separation of clear and cloud conditions for the calibrating period
  calibrating_period$time = as.POSIXct(paste(calibrating_period$Year, "-",
                                             calibrating_period$Month, "-",
                                             calibrating_period$Day, " ",
                                             calibrating_period$Hour, ":",
                                             calibrating_period$Minute, sep=""), tz = timezone)
  top = TOA(latitude_calibrat,calibrating_period$Elev,calibrating_period$time)
  clear_sky_calibrat_DNI = SOLIS(top,calibrating_period$Elev,
                                    aod_380 = rep(threshold_aod_380, length(calibrating_period$Minute)),
                                    aod_500 = rep(threshold_aod_500, length(calibrating_period$Minute)),
                                    w = rep(1, length(calibrating_period$Minute)),
                                    latitude_calibrat)
  clear_data_calibrat = which(calibrating_period$DNI.mod >= clear_sky_calibrat_DNI)
  cloud_data_calibrat = which(calibrating_period$DNI.mod < clear_sky_calibrat_DNI)

# Separate adaptation
  for(data_type in 1:2){

    if(data_type == 1){
      calibrating_period_splitted=calibrating_period[clear_data_calibrat,]
      target_period_splitted=target_period[clear_data_target,]
    }

    if(data_type == 2){
      calibrating_period_splitted=calibrating_period[cloud_data_calibrat,]
      target_period_splitted=target_period[cloud_data_target,]
    }


    # Preprocessing of target period
    if(data_type == 1){print("Preprocessing of clear sky radiation")}
    if(data_type == 2){print("Preprocessing of cloud sky radiation")}
      df_ajustado = pre_process(
      target_period_splitted,latitude_target,z_target,
      calibrating_period_splitted,latitude_calibrat,z_calibrat,
      GHI_threshold,DNI_threshold)

    # Preprocessing of calibrating period
    target_period_splitted$GHI.mod=df_ajustado$GHI_adapted
    target_period_splitted$DNI.mod=df_ajustado$DNI_adapted
    target_period_splitted$DHI.mod=df_ajustado$DHI_adapted

    df_calibrating_ajustado = pre_process(
      calibrating_period_splitted,latitude_calibrat,z_calibrat,
      calibrating_period_splitted,latitude_calibrat,z_calibrat,
      GHI_threshold,DNI_threshold)

    calibrating_period_splitted$GHI.mod = df_calibrating_ajustado$GHI_adapted
    calibrating_period_splitted$DNI.mod = df_calibrating_ajustado$DNI_adapted
    calibrating_period_splitted$DHI.mod = df_calibrating_ajustado$DHI_adapted

    # Site adaptation
    if(data_type == 1){print("Adaptation and postprocessing of clear sky radiation")}
    if(data_type == 2){print("Adaptation and postprocessing  of cloud sky radiation")}
    df_ajustado = adapt_process(
      target_period_splitted,latitude_target,z_target,
      calibrating_period_splitted,latitude_calibrat,z_calibrat,
      GHI_threshold,DNI_threshold)

    if(data_type == 1){adapted_clear_data = df_ajustado}
    if(data_type == 2){adapted_cloudy_data = df_ajustado}

  }

# merge and order adapted data
  adapted_data_above_threshold = data.frame(rbind(adapted_clear_data,adapted_cloudy_data))
  adapted_data_above_threshold = adapted_data_above_threshold[ order(adapted_data_above_threshold$time) , ]

# Data below the solar elevation threshold
  Target_data_below_threshold = Target[which(Target$Elev <= threshold), ]

  df_below_threshold = data.frame(
    Year = Target_data_below_threshold$Year,
    Month = Target_data_below_threshold$Month,
    Day = Target_data_below_threshold$Day,
    Hour = Target_data_below_threshold$Hour,
    Minute = Target_data_below_threshold$Minute,
    GHI_adapted=round(Target_data_below_threshold$GHI.mod,1),
    DNI_adapted=round(Target_data_below_threshold$DNI.mod,1),
    DHI_adapted=round(Target_data_below_threshold$DHI.mod,1),
    Elev = round(Target_data_below_threshold$Elev,2)
  )

  df_above_threshold = data.frame(
    Year = year(adapted_data_above_threshold$time),
    Month = month(adapted_data_above_threshold$time),
    Day = dom(adapted_data_above_threshold$time),
    Hour = hour(adapted_data_above_threshold$time),
    Minute = minute(adapted_data_above_threshold$time),
    GHI_adapted=round(adapted_data_above_threshold$GHI_adapted,1),
    DNI_adapted=round(adapted_data_above_threshold$DNI_adapted,1),
    DHI_adapted=round(adapted_data_above_threshold$DHI_adapted,1),
    Elev = round(adapted_data_above_threshold$Elev,2)
  )

  df_total = rbind(df_above_threshold, df_below_threshold )

  time_stamp = as.POSIXct(paste(df_total$Year, "-",
                                df_total$Month, "-",
                                df_total$Day, " ",
                                df_total$Hour, ":",
                                df_total$Minute, sep=""),
                                tz = timezone)

  df_total=df_total[order(time_stamp),]

  return(df_total)

}
