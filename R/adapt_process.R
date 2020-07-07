#'
#'  Site adaptation of solar irradiance modeled series with coincident ground measurements
#'
#' @param subset_target_period  Dataframe object with solar radiation series to be adapted including time (with same time zone as subset_calibrating_period), the solar irradiance modeled series to be site adapted, along with their clear sky index and solar elevation (degrees)
#' @param latitude_target Site latitude of solar radiation series to be adapted (degrees, +N)
#' @param z_target Site elevation above sea level of solar radiation series to be adapted (m)
#' @param subset_calibrating_period Dataframe object with solar radiation series for calibrating including time (with same time zone as subset_target_period), solar irradiance modeled and measured series, along with modeled clear sky index and solar elevation (degrees)
#' @param latitude_calibrat Site latitude of solar radiation series for calibrating (degrees, +N)
#' @param z_calibrat Site elevation above sea level of solar radiation series for calibrating (m)
#' @param GHI_threshold Upper limit of GHI series (same units that Target). For automatic calulation from observed data, set it to -99
#' @param DNI_threshold Upper limit of DNI series (same units that Target). For automatic calulation from observed data, set it to -99
#'
#' @return Dataframe object including time and site adapted solar irradiance series
#'
#' @export
#' @import glmulti
#' @import solaR
#' @import hyfo




adapt_process <- function(subset_target_period,latitude_target,z_target,
                       subset_calibrating_period,latitude_calibrat,z_calibrat,
                       GHI_threshold,DNI_threshold){

  #--------------------------------------------------------------------------------------------------------------------------
  # Optical air mass
  subset_calibrating_period$air_mass_kasten = m_Kasten(subset_calibrating_period$Elev,z_calibrat)
  subset_target_period$air_mass_kasten = m_Kasten(subset_target_period$Elev,z_target)

  #--------------------------------------------------------------------------------------------------------------------------
  # Top of atmosphere solar irradiance
  subset_calibrating_period$Top_irradiance = TOA(latitude_calibrat,subset_calibrating_period$Elev,subset_calibrating_period$time)
  subset_target_period$Top_irradiance = TOA(latitude_target,subset_target_period$Elev,subset_target_period$time)

  #--------------------------------------------------------------------------------------------------------------------------
  # Non-dimensional parameters
  # Calibrating period
  subset_calibrating_period$KT_model = subset_calibrating_period$GHI.mod / subset_calibrating_period$Top_irradiance
  subset_calibrating_period$KT_measured = subset_calibrating_period$GHI.obs / subset_calibrating_period$Top_irradiance
  subset_calibrating_period$KT_measured_modified = subset_calibrating_period$KT_measured / (1.031*exp( -1.4/(0.9 + 9.4/subset_calibrating_period$air_mass_kasten)) + 0.1)
  subset_calibrating_period$KT_model_modified = subset_calibrating_period$KT_model / (1.031*exp( -1.4/(0.9 + 9.4/subset_calibrating_period$air_mass_kasten)) + 0.1)
  subset_calibrating_period$kd_measured =  subset_calibrating_period$DHI.obs /subset_calibrating_period$GHI.obs
  subset_calibrating_period$kd_model = subset_calibrating_period$DHI.mod/subset_calibrating_period$GHI.mod
  subset_calibrating_period = subset_calibrating_period[which(subset_calibrating_period$kd_measured > 0),]
  subset_calibrating_period = subset_calibrating_period[which(subset_calibrating_period$kd_measured <= 1),]

  # Target period
  subset_target_period$KT_model = subset_target_period$GHI.mod / subset_target_period$Top_irradiance
  subset_target_period$kd_model = subset_target_period$DHI.mod /subset_target_period$GHI.mod
  subset_target_period$KT_model_modified = subset_target_period$KT_model / (1.031*exp( -1.4/(0.9 + 9.4/subset_target_period$air_mass_kasten)      )   + 0.1)

  #--------------------------------------------------------------------------------------------------------------------------
  # The local adaptation only applies to solar Elevs above 1 degree:
  threshold = 1
  subset_target_period_low_elev = subset_target_period[which(subset_target_period$Elev <= threshold),]
  subset_calibrating_period_high_elev = subset_calibrating_period[which(subset_calibrating_period$Elev > threshold),]
  subset_target_period_high_elev = subset_target_period[which(subset_target_period$Elev > threshold),]

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$                                GHI                                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  frc_ajuste = subset_target_period_high_elev[, c(which(names(subset_target_period_high_elev) == "time"),
                                           which(names(subset_target_period_high_elev) == "GHI.mod"))]

  frc_ajuste$time = as.Date(frc_ajuste$time)
  hindcast_ajuste = subset_calibrating_period_high_elev[, c(which(names(subset_calibrating_period_high_elev) == "time"),
                                                     which(names(subset_calibrating_period_high_elev) == "GHI.mod"))]
  hindcast_ajuste$time = as.Date(hindcast_ajuste$time)
  obs_ajuste = subset_calibrating_period_high_elev[, c(which(names(subset_calibrating_period_high_elev) == "time"),
                                                which(names(subset_calibrating_period_high_elev) == "GHI.obs"))]
  obs_ajuste$time = as.Date(obs_ajuste$time)
  GHI_adaptation_output <- biasCorrect(frc_ajuste,
                                           hindcast_ajuste,
                                           obs_ajuste,
                                           method = 'eqm',
                                           extrapolate = "no")
  GHI_Locally_adapted = GHI_adaptation_output$GHI.mod

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$                                DNI                                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  frc_ajuste = subset_target_period_high_elev[, c(which(names(subset_target_period_high_elev) == "time"),
                                           which(names(subset_target_period_high_elev) == "DNI.mod"))]
  frc_ajuste$time = as.Date(frc_ajuste$time)
  hindcast_ajuste = subset_calibrating_period_high_elev[, c(which(names(subset_calibrating_period_high_elev) == "time"),
                                                     which(names(subset_calibrating_period_high_elev) == "DNI.mod"))]
  hindcast_ajuste$time = as.Date(hindcast_ajuste$time)
  obs_ajuste = subset_calibrating_period_high_elev[, c(which(names(subset_calibrating_period_high_elev) == "time"),
                                                which(names(subset_calibrating_period_high_elev) == "DNI.obs"))]
  obs_ajuste$time = as.Date(obs_ajuste$time)
  DNI_adaptation_output <- biasCorrect(frc_ajuste,
                                               hindcast_ajuste,
                                               obs_ajuste,
                                               method = 'eqm',
                                               #method = 'scaling',
                                               extrapolate = "no")
  DNI_Fit = DNI_adaptation_output$DNI.mod

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$                                DNI                                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  frc_ajuste = subset_target_period_high_elev[, c(which(names(subset_target_period_high_elev) == "time"),
                                           which(names(subset_target_period_high_elev) == "DHI.mod"))]
  frc_ajuste$time = as.Date(frc_ajuste$time)
  hindcast_ajuste = subset_calibrating_period_high_elev[, c(which(names(subset_calibrating_period_high_elev) == "time"),
                                                     which(names(subset_calibrating_period_high_elev) == "DHI.mod"))]
  hindcast_ajuste$time = as.Date(hindcast_ajuste$time)
  obs_ajuste = subset_calibrating_period_high_elev[, c(which(names(subset_calibrating_period_high_elev) == "time"),
                                                which(names(subset_calibrating_period_high_elev) == "DHI.obs"))]
  obs_ajuste$time = as.Date(obs_ajuste$time)
  DHI_adaptation_output <- biasCorrect(frc_ajuste,
                                               hindcast_ajuste,
                                               obs_ajuste,
                                               method = 'eqm',
                                               #method = 'scaling',
                                               extrapolate = "no") # este es el mejor
  DHI_Fit = DHI_adaptation_output$DHI.mod

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$                Adapted series                                 $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  df_daytime <- data.frame(
    time = subset_target_period_high_elev$time,
    GHI_adapted = as.numeric(GHI_Locally_adapted),
    DNI_adapted = as.numeric(DNI_Fit),
    DHI_adapted = as.numeric(DHI_Fit),
    Elev = subset_target_period_high_elev$Elev
  )

  #################################################################################################################
  #  Postprocessing

  df_daytime_original = df_daytime
  if(GHI_threshold == -99){GHI_threshold = 1.05*max(subset_calibrating_period$GHI.obs)}
  if(DNI_threshold == -99){DNI_threshold = 1.05*max(subset_calibrating_period$DNI.obs)}

  df_daytime = post_process(df_daytime_original,
                            subset_target_period_high_elev,
                            GHI_threshold, DNI_threshold)


  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$          Merge of adapted and low elevation dataframes       $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  df_low_elev<- data.frame(
    time = subset_target_period_low_elev$time,
    GHI_adapted = subset_target_period_low_elev$GHI.mod,
    DNI_adapted = subset_target_period_low_elev$DNI.mod,
    DHI_adapted = subset_target_period_low_elev$DHI.mod,
    Elev = subset_target_period_low_elev$Elev
  )

  sal = df_daytime
  if(length(df_low_elev$time) > 0){sal = rbind(df_daytime,df_low_elev)}
  sal = sal[order(year(sal$time)*1000 + doy(sal$time) + (round(hour(sal$time)/24, 2)) ),]

  return(sal)

}
