#'
#'  Postprocessing of adapted solar irradiance
#'
#' @param subset_target_period_high_elev  Dataframe object with daytime modeled solar radiation series including time (with same time zone as df_daytime)
#' @param df_daytime Dataframe object with daytime adapted solar radiation series including time (with same time zone as subset_target_period_high_elev)
#' @param GHI_threshold GHI threshold value, in the same units that modeled and adapted datasets. Default value is -99
#' @param DNI_threshold GHI threshold value, in the same units that modeled and adapted datasets. Default value is -99
#'
#' @return Dataframe object including time and site adapted solar irradiance series without inconsistencies
#'
#' @export


post_process <- function(df_daytime,subset_target_period_high_elev, GHI_threshold, DNI_threshold){

# Correction 1: negative solar irradiance values
  # Negative GHI values
  if(length(which(df_daytime$GHI_adapted < 0)) > 0){
    negative_values = which(df_daytime$GHI_adapted < 0)
    df_daytime$GHI_adapted[negative_values] = subset_target_period_high_elev$GHI.mod[negative_values]
    df_daytime$DNI_adapted[negative_values] = subset_target_period_high_elev$DNI.mod[negative_values]
    df_daytime$DHI_adapted[negative_values] = subset_target_period_high_elev$DHI.mod[negative_values]
  }
  # Negative DNI values
  if(length(which(df_daytime$DNI_adapted < 0)) > 0){
    negative_values = which(df_daytime$DNI_adapted < 0)
    df_daytime$GHI_adapted[negative_values] = subset_target_period_high_elev$GHI.mod[negative_values]
    df_daytime$DNI_adapted[negative_values] = subset_target_period_high_elev$DNI.mod[negative_values]
    df_daytime$DHI_adapted[negative_values] = subset_target_period_high_elev$DHI.mod[negative_values]
  }
  # Negative DHI values
  if(length(which(df_daytime$DHI_adapted < 0)) > 0){
    negative_values = which(df_daytime$DHI_adapted < 0)
    df_daytime$GHI_adapted[negative_values] = subset_target_period_high_elev$GHI.mod[negative_values]
    df_daytime$DNI_adapted[negative_values] = subset_target_period_high_elev$DNI.mod[negative_values]
    df_daytime$DHI_adapted[negative_values] = subset_target_period_high_elev$DHI.mod[negative_values]
  }

# Correction 2: Extremely high solar irradiance values
  # Extremely high GHI values
  if(length(which(df_daytime$GHI_adapted > GHI_threshold)) > 0){
    extremely_high = which(df_daytime$GHI_adapted > GHI_threshold)
    df_daytime$GHI_adapted[extremely_high] = subset_target_period_high_elev$GHI.mod[extremely_high]
    df_daytime$DNI_adapted[extremely_high] = subset_target_period_high_elev$DNI.mod[extremely_high]
    df_daytime$DHI_adapted[extremely_high] = subset_target_period_high_elev$DHI.mod[extremely_high]
  }

  # Extremely high DNI values
  if(length(which(df_daytime$DNI_adapted > DNI_threshold)) > 0){
    extremely_high = which(df_daytime$DNI_adapted > DNI_threshold)
  df_daytime$GHI_adapted[extremely_high] = subset_target_period_high_elev$GHI.mod[extremely_high]
  df_daytime$DNI_adapted[extremely_high] = subset_target_period_high_elev$DNI.mod[extremely_high]
  df_daytime$DHI_adapted[extremely_high] = subset_target_period_high_elev$DHI.mod[extremely_high]
  }

  # Extremely high DHI values
  if(length(which(df_daytime$DHI_adapted > 700))){
    extremely_high = which(df_daytime$DHI_adapted > 700)
    df_daytime$GHI_adapted[extremely_high] = subset_target_period_high_elev$GHI.mod[extremely_high]
    df_daytime$DNI_adapted[extremely_high] = subset_target_period_high_elev$DNI.mod[extremely_high]
    df_daytime$DHI_adapted[extremely_high] = subset_target_period_high_elev$DHI.mod[extremely_high]
  }


# Correction 3: Inconsistent data
  # Solar zenith angle < 75
  inconsistent_1 = intersect(
    which(df_daytime$Elev > 15),
    which(df_daytime$GHI_adapted /
            (df_daytime$DHI_adapted +
               df_daytime$DNI_adapted*sin(df_daytime$Elev*pi/180)) > 1.08 |
            df_daytime$GHI_adapted /
            (df_daytime$DHI_adapted +
               df_daytime$DNI_adapted*sin(df_daytime$Elev*pi/180)) < 0.92))

  if(length(inconsistent_1) > 0){
    df_daytime$GHI_adapted[inconsistent_1] = subset_target_period_high_elev$GHI.mod[inconsistent_1]
    df_daytime$DNI_adapted[inconsistent_1] = subset_target_period_high_elev$DNI.mod[inconsistent_1]
    df_daytime$DHI_adapted[inconsistent_1] = subset_target_period_high_elev$DHI.mod[inconsistent_1]
  }

  inconsistent_2 = intersect(
    which(df_daytime$Elev <= 15),
    which(df_daytime$GHI_adapted /
            (df_daytime$DHI_adapted +
               df_daytime$DNI_adapted*sin(df_daytime$Elev*pi/180)) > 1.15 |
            df_daytime$GHI_adapted /
            (df_daytime$DHI_adapted +
               df_daytime$DNI_adapted*sin(df_daytime$Elev*pi/180)) < 0.85))

  if(length(inconsistent_2) > 0){
    df_daytime$GHI_adapted[inconsistent_2] = subset_target_period_high_elev$GHI.mod[inconsistent_2]
    df_daytime$DNI_adapted[inconsistent_2] = subset_target_period_high_elev$DNI.mod[inconsistent_2]
    df_daytime$DHI_adapted[inconsistent_2] = subset_target_period_high_elev$DHI.mod[inconsistent_2]
  }


# Correction 4: Suspicious data
  # correccion_errores
  suspicius = intersect(which(subset_target_period_high_elev$GHI.mod > 200),
                       which(df_daytime$GHI_adapted < 10))
  if(length(suspicius) > 0){
    df_daytime$GHI_adapted[suspicius] = subset_target_period_high_elev$GHI.mod[suspicius]
    df_daytime$DNI_adapted[suspicius] = subset_target_period_high_elev$DNI.mod[suspicius]
    df_daytime$DHI_adapted[suspicius] = subset_target_period_high_elev$DHI.mod[suspicius]
  }

  return(df_daytime)

}
