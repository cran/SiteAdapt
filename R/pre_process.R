#'
#'  Preprocessing of solar irradiance series Site adaptation
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


pre_process <- function(subset_target_period,
                        latitude_target,
                        z_target,
                        subset_calibrating_period,
                        latitude_calibrat,
                        z_calibrat,
                        GHI_threshold,
                        DNI_threshold){


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
  subset_calibrating_period$KT_measured_modified = subset_calibrating_period$KT_measured / (1.031*exp( -1.4/(0.9 + 9.4/subset_calibrating_period$air_mass_kasten)      )   + 0.1)
  subset_calibrating_period$KT_model_modified = subset_calibrating_period$KT_model / (1.031*exp( -1.4/(0.9 + 9.4/subset_calibrating_period$air_mass_kasten)      )   + 0.1)
  subset_calibrating_period$kd_measured = (subset_calibrating_period$GHI.obs - subset_calibrating_period$DNI.obs * sin(subset_calibrating_period$Elev*pi/180))/subset_calibrating_period$GHI.obs
  subset_calibrating_period$kd_model = (subset_calibrating_period$GHI.mod - subset_calibrating_period$DNI.mod * sin(subset_calibrating_period$Elev*pi/180))/subset_calibrating_period$GHI.mod
  subset_calibrating_period = subset_calibrating_period[which(subset_calibrating_period$kd_measured > 0),]
  subset_calibrating_period = subset_calibrating_period[which(subset_calibrating_period$kd_measured <= 1),]

  # Target period
  subset_target_period$KT_model = subset_target_period$GHI.mod / subset_target_period$Top_irradiance
  subset_target_period$kd_model = (subset_target_period$GHI.mod - subset_target_period$DNI.mod * sin(subset_target_period$Elev*pi/180))/subset_target_period$GHI.mod
  subset_target_period$KT_model_modified = subset_target_period$KT_model / (1.031*exp( -1.4/(0.9 + 9.4/subset_target_period$air_mass_kasten)      )   + 0.1)

  # DHI
  subset_target_period$DHI.mod = (subset_target_period$GHI.mod - subset_target_period$DNI.mod * sin(subset_target_period$Elev*pi/180))

  if(length(which(subset_target_period$DHI.mod < 0)) > 0){
    DHI_modelo_inconsistente = which(subset_target_period$DHI.mod < 0)
    subset_target_period$DHI.mod[DHI_modelo_inconsistente] = subset_target_period$GHI.mod[DHI_modelo_inconsistente]
  }

  #--------------------------------------------------------------------------------------------------------------------------
  # The local adaptation only applies to solar Elevs above 1 degree:
  subset_target_period_low_elev = subset_target_period[which(subset_target_period$Elev <= 1),]
  subset_calibrating_period_high_elev = subset_calibrating_period[which(subset_calibrating_period$Elev > 1),]
  subset_target_period_high_elev = subset_target_period[which(subset_target_period$Elev > 1),]

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$                                GHI and DHI                          $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  GHI_adaptation_proc <-  glm(KT_measured ~ KT_model + kd_model + kc + Elev,data=subset_calibrating_period_high_elev)
  Best_GHI_adaptation <- glmulti(GHI_adaptation_proc,level = 2,crit="AIC", family = gaussian(link = "identity"), confsetsize=100 ,plotty = F, report = F)
  model_predicting_kt <- Best_GHI_adaptation@objects[[1]];vars <- names(coef(model_predicting_kt))[-1]
  kt_Locally_adapted = as.numeric(predict(model_predicting_kt, subset_target_period_high_elev, family = gaussian(link = "identity")))
  GHI_Locally_adapted = kt_Locally_adapted*subset_target_period_high_elev$Top_irradiance
  GHI_Locally_adapted = kt_Locally_adapted*subset_target_period_high_elev$Top_irradiance

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$                                DNI                                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  adapted_kt_for_subset_calibrating_period <- predict(model_predicting_kt, subset_calibrating_period_high_elev, family = gaussian(link = "identity"))
  subset_calibrating_period_high_elev$non_dim_kt = as.numeric(adapted_kt_for_subset_calibrating_period)/(1.031*exp( -1.4/(0.9 + 9.4/subset_calibrating_period_high_elev$air_mass_kasten))+ 0.1)
  subset_target_period_high_elev$non_dim_kt = kt_Locally_adapted/(1.031*exp( -1.4/(0.9 + 9.4/subset_target_period_high_elev$air_mass_kasten))+ 0.1)
  DNI_adaptation_proc <-  glm(kd_measured ~ poly(non_dim_kt, 4) +  Elev + kd_model + kc, data=subset_calibrating_period_high_elev)
  Best_DNI_adaptation <- glmulti(DNI_adaptation_proc,level = 2, crit="AIC", family =  "gaussian", confsetsize=100, plotty = FALSE, report = FALSE)
  model <- Best_DNI_adaptation@objects[[1]];vars <- names(coef(model))[-1]
  salida_kd <- predict(model, subset_target_period_high_elev)
  DHI_Fit = as.numeric(salida_kd)*GHI_Locally_adapted
  DNI_Fit = (GHI_Locally_adapted - DHI_Fit) / sin(subset_target_period_high_elev$Elev*pi/180)

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$                Adapted series                                 $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  df_daytime <- data.frame(
    Time = subset_target_period_high_elev$time,
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
  df_daytime = post_process(df_daytime_original,subset_target_period_high_elev,GHI_threshold, DNI_threshold)


  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$          Merge of adapted and low elevation dataframes       $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  df_low_elev<- data.frame(
    Time = subset_target_period_low_elev$time,
    GHI_adapted = subset_target_period_low_elev$GHI.mod,
    DNI_adapted = subset_target_period_low_elev$DNI.mod,
    DHI_adapted = subset_target_period_low_elev$DHI.mod
    )

  sal = df_daytime
  if(length(df_low_elev$Time) > 0){sal = rbind(df_daytime,df_low_elev)}
  sal = sal[order(year(sal$Time)*1000 + doy(sal$Time) + (round(hour(sal$Time)/24, 2)) ),]

  return(sal)

}
