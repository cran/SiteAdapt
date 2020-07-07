#'
#'  DNI under clear sky conditions (Solar Energy, 82(8), 758-762)
#'
#' @param top  Solar irradiance at the top of atosphere
#' @param Sun_elev Sun elevation angle (degrees) above horizon
#' @param aod_380 Aerosol optical depth at 380 nm (dimensionless)
#' @param aod_500 Aerosol optical depth at 500 nm (dimensionless)
#' @param w Atmospheric water vapor content (cm)
#' @param site_elevation Site elevation above sea level (m)
#'
#' @return Dataframe object including time and site adapted solar irradiance series
#'
#' @export


SOLIS <- function(top, Sun_elev, aod_380, aod_500, w, site_elevation){

# Kasten optic air mass
gamma = Sun_elev*pi/180
Kasten_Opt_Air_mass = m_Kasten(gamma,site_elevation)

# Parameters
Ibn=0
aod_700 = 0.27583*aod_380 + 0.35*aod_500
I00 = 1.08*(w)^(0.0051)
I01 = 0.97*(w)^(0.032)
I02 = 0.12*(w)^(0.56)

I_prima_0 = top * (I02 * (aod_700)^2 + I01 * aod_700 + I00 + 0.071*(-site_elevation/8435.2))

# tau_b
tb1 = 1.82+0.056*log(w) + 0.0071*((log(w))^2)
tb0 = 0.33+0.045*log(w) + 0.0096*((log(w))^2)
tbp = 0.0089*w + 0.13
tau_b = tb1 * aod_700 + tb0 + tbp*(-site_elevation/8435.2)

# b
b1 = 0.00925*((aod_700)^2) + 0.0148*aod_700 - 0.0172
b0 = -0.7565*((aod_700)^2) + 0.5057*aod_700 + 0.4557
b = b1*log(w) + b0

Ibn = I_prima_0 * exp(- ( tau_b / (sin(gamma)^b  )))

return(Ibn)

}


