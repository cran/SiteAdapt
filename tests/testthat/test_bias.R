context("Bias of adapted solar radiation")


# A site located in the the Namib Desert of Namibia (Gobabeb, GOB) is selected

# - latitude:   23.5614 S
# - Longitude:  15.0420 E
# - Elevation:  407.0 m asl

# Load calibration and modeled datasets
data(calibration_2016) # Measured from BSRN
data(target_2013_2016) # Provided by CAMS-RAD service

observed_2013_2016$time = as.POSIXct(
  paste(observed_2013_2016$Year, "-",
        observed_2013_2016$Month, "-",
        observed_2013_2016$Day, " ",
        observed_2013_2016$Hour, ":",
        observed_2013_2016$Minute, sep=""),
  tz ="UTC")

target_2013_2016$time = as.POSIXct(
  paste(target_2013_2016$Year, "-",
        target_2013_2016$Month, "-",
        target_2013_2016$Day, " ",
        target_2013_2016$Hour, ":",
        target_2013_2016$Minute, sep=""),
  tz ="UTC")

# Apply the site adaptation procedure
site_adapted_series = site_adapt(
  Target = target_2013_2016,
  latitude_target = -23.5614, # Latitude of target site
  z_target = 407.0, # Elevation of target site
  Calibration = calibration_2016,
  latitude_calibrat = -23.5614, # Same location of target period
  z_calibrat = 407.0, # Same location of target period
  timezone = "UTC",
  GHI_threshold = -99, # The threshold is calculated from observed data
  DNI_threshold = -99) # The threshold is calculated from observed data


# Load measured data for evaluating the site adaptation performance:
data(observed_2013_2016)

# Merge datasets

site_adapted_series$time = as.POSIXct(
  paste(site_adapted_series$Year, "-",
        site_adapted_series$Month, "-",
        site_adapted_series$Day, " ",
        site_adapted_series$Hour, ":",
        site_adapted_series$Minute, sep=""),
  tz ="UTC")

meas_model = merge(observed_2013_2016[,6:9],
                   target_2013_2016[,c(6:9,11)],
                   by = "time", all = FALSE )

meas_model_adapt = merge(meas_model,
                         site_adapted_series[,6:10],
                         by = "time", all = FALSE )

# Statistical indicators
library(hydroGOF)

bias_DNI_adapted = pbias(meas_model_adapt$DNI_adapted,meas_model_adapt$DNI.obs)
bias_GHI_adapted = pbias(meas_model_adapt$GHI_adapted,meas_model_adapt$GHI.obs)

test_that("Adapted DNI bias is lower than 1%", {
  expect_that(floor(abs(bias_DNI_adapted)), equals(0))
})

test_that("Adapted GHI bias is lower than 1%", {
  expect_that(floor(abs(bias_GHI_adapted)), equals(0))
})

