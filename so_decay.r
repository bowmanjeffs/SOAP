setwd('~/ducklow_lab/palmer_2015')

## This script is based on the calibration portion of the regular SO analysis
## which is reflected in the nomenclature.  Unlike the regular calibration however,
## the parameters for each decay experiment should be written into a tsv file in the order:

## data file name
## time standard is read on spec
## time tube is placed in sample
## first reading at 240
## final reading at 240
## aliquot size in ml

### file names ####

data.path <- 'data/so_data/' # path to the parameter and data files
calib.file <- '15.12.51_decay.txt' # file name for calibration results output
param.file <- '15.12.21_decay_data.txt' # name of file with experimental parameters

#### calibration ####

### enter calibration parameters here ###

flow.rate <- 7.5 # flow rate in ml min-1

calib.params <- read.table(paste0(data.path, param.file), header = T, stringsAsFactors = F)

### set up matrix to hold decay experiment results ###

calib.results <- matrix(ncol = 11, nrow = length(row.names(calib.params)))
row.names(calib.results) <- c(calib.params$decay.file)
colnames(calib.results) <- c('L0', 'T0', 'R2', 'k', 'half.life', 'start.buffer', 'end.buffer', 'start.std', 'end.std', 'start.sod', 'end.sod')

for(c in 1:length(row.names(calib.params))){
  
  ### manually entered parameters ###
  
  spec.time <- as.numeric(calib.params[c, 2]) # time standard was read on spec
  sample.time <- as.numeric(calib.params[c, 3]) # time ko2 goes on the instrument
  calib.240 <- as.numeric(calib.params[c, 4]) - as.numeric(calib.params[c, 5]) # abs at 240
  aliquot <- as.numeric(calib.params[c, 6]) # size of standard aliquot in L
  
  ## read in times
  
  calib <- read.table(paste0(data.path, calib.params[c, 1]), skip = 1)
  calib$V1 <- calib$V1 * 2 # convert from events to seconds
  
  colnames(calib) <- c('s', 'lum')
  
  ### generate plot and pick points ###
  
  plot(calib$lum ~ calib$s,
       type = 'l')
  
  ## pick points along stable buffer line
  print('identify two stable buffer points')
  buffer.start <- calib$s[identify(calib$s, calib$lum, n = 1, plot = T)]
  buffer.end <- calib$s[identify(calib$s, calib$lum, n = 1, plot = T)]
  
  ## pick peak of SO spike and point along decay curve
  print('identify peak of SO spike and point along decay curve')
  ko2_start <- calib$s[identify(calib$s, calib$lum, n = 1, plot = T)]
  ko2_end <- calib$s[identify(calib$s, calib$lum, n = 1, plot = T)]
  
  ## pick two points along post-SOD baseline
  print('identify two points along post-SOD baseline')
  sod.start <- calib$s[identify(calib$s, calib$lum, n = 1, plot = T)]
  sod.end <- calib$s[identify(calib$s, calib$lum, n = 1, plot = T)] 
  
  ## average values for times
  
  buffer.mean <- mean(calib$lum[calib$s >= buffer.start & calib$s <= buffer.end])
  buffer.sd <- sd(calib$lum[calib$s >= buffer.start & calib$s <= buffer.end])
  
  sod.mean <- mean(calib$lum[calib$s >= sod.start & calib$s <= sod.end])
  sod.sd <- sd(calib$lum[calib$s >= sod.start & calib$s <= sod.end])
  
  ### calculated parameters ###
  
  ext.time <- sample.time - spec.time # delay between spec reading and putting ko2 on instrument
  calib.time <- ko2_start - ext.time # time between ko2 getting on instrument and ko2 reaching flow cell
  baseline.drift <- buffer.mean - sod.mean
  ko2.stock <- calib.240 / 2183 # 2183 = molar absorptivity at 240 Bielski et al 1985
  std.conc <- ko2.stock * aliquot / (0.01 + aliquot) * 1e12 # concentration of standard at T0 in pMol
  
  ## linearized decay curve calculation
  
  log.decay <- calib[calib$s >= ko2_start & calib$s <= ko2_end,]
  log.decay$s <- log.decay$s - calib.time
  log.decay$lum <- log.decay$lum - buffer.mean
  linear.decay <- lm(log(log.decay$lum) ~ log.decay$s)
  linear.decay.sm <- summary(linear.decay)
  
  plot(log(log.decay$lum) ~ log.decay$s,
       xlim = c(0, 150),
       ylim = c(7, 15))
  abline(linear.decay)
  
  l_0 <- exp(linear.decay$coefficients[1])
  k <- linear.decay$coefficients[2] * -1
  half.life <- log(2) / k / 60 # half life is in minutes
  
  ### output calibration results ###
  
  cal.slope <- std.conc / l_0
  calib.results[c,] <- c(l_0, std.conc, linear.decay.sm$r.squared, k, half.life, buffer.start, buffer.end, ko2_start, ko2_end, sod.start, sod.end)
  
}

write.table(calib.results, paste0(data.path, calib.file), quote = F, row.names = T)
