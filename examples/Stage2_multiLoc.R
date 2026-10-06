################################################################################
#
# Second stage: multi-location studies
#
# Example with Italian urban data
#
################################################################################

# Analysis
library(dlnm); library(splines); library(mixmeta)

# Convenience
library(tidyverse); library(doParallel); library(foreach)

#-------------------------
# Parameters of the analysis
#-------------------------

#----- DLNM

# Exposure dimension
varfun <- "bs"
varper <- c(10,75,90)
vardegree <- 2

# Lag dimension
maxlag <- 21
lagfun <- "ns"
lagknots <- logknots(maxlag, 3)

#----- For plotting and other

# Percentiles to evaluate
tper <- c(seq(0, 1, by = .1), 2:98, seq(99, 100, by = .1))

#-------------------------
# Prepare data
#-------------------------

# Load full data
italy <- read.csv("data/italy_ts.csv", colClasses = c("date" = "Date"))

# Split by city
dlist <- split(italy[,-1], italy$city_code)

# Get decription of cities
metadata <- read.csv("data/italy_meta.csv")

#-------------------------
# Run first-stage across cities
#-------------------------

#----- Prepare parallelisation
ncores <- detectCores()
cl <- makeCluster(max(1, ncores - 2))
registerDoParallel(cl)

#----- Perform first stage on each city

# Loop over cities, splitting the data table
stage1res <- foreach(dat = dlist,
  .packages = c("dlnm", "splines", "foreach", "mixmeta", "dplyr")) %dopar% 
{

  # Define crossbasis
  argvar <- list(fun = varfun, degree = vardegree, 
    knots = quantile(dat$tmean, varper / 100, na.rm = T))
  cb <- crossbasis(dat$tmean, lag = maxlag, argvar = argvar,
    arglag = list(fun = lagfun, knots = lagknots))
    
  # Run model
  res <- glm(deaths ~ cb + dow + ns(date, df = 7 * length(unique(year))), 
    dat, family = quasipoisson)
    
  # Reduce coefficients to overall cumulative
  # This is to meta-analyse a lower number of coefficients
  redall <- crossreduce(cb, res, cen = median(dat$tmean, na.rm = T))
    
  # Output
  list(coefs = coef(redall), vcov = vcov(redall), 
    tdist = quantile(dat$tmean, tper / 100))
}

# Stop parallel
stopCluster(cl)

#-------------------------
# Second stage
#-------------------------

# Extract first-stage results
coefs <- sapply(stage1res, "[[", "coefs") |> t()
vcovs <- lapply(stage1res, "[[", "vcov")
cities <- names(stage1res)

#----- Fixed effect meta-analysis

# Fit model
fixmod <- mixmeta(coefs, S = vcovs, method = "fixed")

# Extract overall effect
fixeff <- predict(fixmod, vcov = T)[[1]]

#----- Random effect meta-analysis

# Fit model
ranmod <- mixmeta(coefs, S = vcovs, random = ~ 1|city_code, data = metadata)

# Extract overall effect
raneff <- predict(ranmod, vcov = T)[[1]]

# Extract blups
ranblups <- blup(ranmod, vcov = T)

#----- Meta-regression

# Fit random-effect meta-regression
metamod <- mixmeta(coefs ~ gdp, S = vcovs, random = ~ 1|city_code, 
  data = metadata)

# Extract predictions at two values
preds <- predict(metamod, 
  data.frame(gdp = quantile(metadata$gdp, c(.10, .90))),
  vcov = T)

# Extract improved estimates (BLUPs)
blup_coefs <- blup(metamod, vcov = T)

#-------------------------
# Save results
#-------------------------

save.image("data/stage2.RData")