#'
#'  Calculation of relative air mass based on Kasten parametrization
#'
#' @param Sun_elev Sun elevation angle (degrees) above horizon
#' @param z Site elevation above sea level (m)
#'
#' @return Relative air mass based on Kasten parametrization
#'
#' @export


m_Kasten <- function(Sun_elev,z){

ier=0;gamma=Sun_elev*pi/180
# sun elevation correction for refraction
c1 <- 0.061359;c2 <- 0.1594;c3 <- 1.1230;c4 <- 0.065656;c5 <- 28.9344;c6 <- 277.3971
gamma_corr <- c1*(c2+c3*gamma+c4*gamma*gamma)/(1.+c5*gamma+c6*gamma*gamma);gamma <- gamma+gamma_corr
gamma_deg <- Sun_elev;
a <- 0.50572
b <- 6.07995 #degrees
c <- -1.6364
m = exp(-z/8434.5) /(sin(gamma)+a*(Sun_elev+b)^c)

return(m)
}
