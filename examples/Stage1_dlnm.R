################################################################################
#
# First stage: distributed-lag (nonlinear) models
#
# Example with London data
#
################################################################################

#----------------------
# Useful packages
#----------------------

# Analysis
library(dlnm); library(splines)

# Convenience
library(tidyverse)

#----------------------
# Some parameters
#----------------------

# ns degrees of freedom for each year
dfseas <- 7

# Maximum lag for PM10
maxlag <- 7

#----------------------
# Prepare data
#----------------------

# Read data
london <- read.csv("data/london.csv") |>
  mutate(date = as.Date(date, format = "%d/%m/%Y")) |>
  subset(year(date) < 2002)

# Add day-of-week variabe
london <- mutate(london, dow = factor(wday(date)))

# Get the number of years and df for spline
ny <- length(unique(year(london$date)))

#----------------------
# First Poisson regression
# Same day PM10
#----------------------

# Fit model
m0 <- glm(death ~ pm10 + ns(london$date, df = ny * dfseas) + dow, data = london, 
  family = "quasipoisson", na.action = "na.exclude")

# Get coefficients and confidence intervals
coef0 <- coef(m0)
ci0 <- confint(m0)

# Get the time effect
preddf <- mutate(london, pm10 = mean(pm10, na.rm = T), dow = dow[1])
timeff <- predict(m0, type = "response", newdata = preddf)

#----------------------
# Lagged association
#----------------------

# Define lagged exposure
xlags <- sapply(0:maxlag, lag, x = london$pm10)
colnames(xlags) <- sprintf("Lag%i", 0:maxlag)

# Fit model
m1 <- glm(death ~ xlags + ns(date, ny * 7) + dow, data = london,
  family = "quasipoisson", na.action = "na.exclude")

# Extract coefs and confidence intervals
coef1 <- coef(m1)
ci1 <- confint(m1)

#----------------------
# Distributed lag model
#----------------------

# We use the package dlnm to create the basis
xdl <- crossbasis(london$pm10, lag = maxlag, 
  arglag = list(fun = "ns", knots = logknots(maxlag, 1)))

# Fit model
m2 <- glm(death ~ xdl + ns(date, ny * 7) + dow, data = london,
  family = "quasipoisson", na.action = "na.exclude")

# Extract lag-response function
lrf2 <- crosspred(xdl, m2, at = 1, cen = 0)

#----------------------
# Distributed lag nonlinear model
# On temperature
#----------------------

# Now create crossbasis by also adding a nonlinear spec on var dimensions
xcb <- crossbasis(london$tmean, lag = maxlag, 
  argvar = list(fun = "bs", df = 4),
  arglag = list(fun = "bs", knots = logknots(maxlag, 2)))

# Fit model
m3 <- glm(death ~ xcb + ns(date, ny * 7) + dow, data = london,
  family = "quasipoisson", na.action = "na.exclude")

# Extract full surface
surf3 <- crosspred(xcb, m3)

#----------------------
# Save results
#----------------------

save.image("data/stage1.RData")