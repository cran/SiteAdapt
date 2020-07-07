
#'
#'  Calculation of the Top of Atmosphere (TOA) solar irradiance on a horizontal plane
#'
#' @param latitude Site latitude (degrees, +N)
#' @param Sun_elev Sun elevation angle (degrees) above horizon
#' @param Time_Stamp Time series (object of class "POSIXct")
#'
#' @return Top of Atmosphere (TOA) Solar irradiance on a horizontal plane
#'
#' @export
#' @import solaR

TOA <- function(latitude,Sun_elev,Time_Stamp){

  solDs<-fSolD(latitude,Time_Stamp)
  yyyy_mm_dd_df=data.frame(yyyy_mm_dd=year(Time_Stamp)*1000 + doy(Time_Stamp))
  yyyy_mm_dd_unico=unique(yyyy_mm_dd_df$yyyy_mm_dd)
  Eo=as.numeric(solDs[,2])
  TOA = Eo*1367
  df_Eo=data.frame(Eo=Eo, yyyy_mm_dd=yyyy_mm_dd_unico, TOA = TOA)
  df_todo_Eo=merge(df_Eo,yyyy_mm_dd_df, by="yyyy_mm_dd",all=TRUE )
  Top = df_todo_Eo$TOA*sin(Sun_elev*pi/180)

return(Top)
}
