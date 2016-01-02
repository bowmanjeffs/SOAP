setwd('~/ducklow_lab/palmer_2015')

## This script is for analyzing runs with multiple experiments in each run.  For example, if you ran
## light and dark treatments sequentially in single run use this.  Currently set up for 4 experiments
## per run, easily modified for fewer or more.

### file names ####

data.path <- 'data/so_data/' # path to the parameter and data files
assay.file <- '15.12.16_so_assays.txt' # file where data will be written to
calib.file <- '15.12.16_calibration.txt' # file name for calibration results output
data.files <- read.table('15.12.16_data_files.txt', header = F, stringsAsFactors = F) # single column list of data files to be analyzed

cal.file.1 <- '15.12.16_cal1.txt' # first calibration file
cal.file.2 <- '15.12.16_cal2.txt' # second calibration file
cal.file.3 <- '15.12.16_cal3.txt' # third calibration file

#### experiment names ####

sample1 <- 'dark.0.2'
sample2 <- 'dark'
sample3 <- 'light.0.2'
sample4 <- 'light'

#### calibration ####

### enter calibration parameters here ###

flow.rate <- 6.5 # flow rate in ml min-1

## set up matrix to hold calibration parameters

calib.n <- 3 # number of calibration points
calib.params <- matrix(nrow = calib.n, ncol = 6)
colnames(calib.params) <- c('calib.file', 'spec.time', 'sample.time', 'calib.240.1', 'calib.240.2', 'aliquot')

## parameters should be entered in order:

# relative path to data file
# time standard is read on spec
# time tube is placed in sample
# first reading at 240
# final reading at 240
# aliquot size in ml

calib.params [1,] <- c(paste0(data.path, cal.file.1),
                       74,
                       109,
                       0.09,
                       0.057,
                       0.000005)
calib.params [2,] <- c(paste0(data.path, cal.file.2),
                       99,
                       137,
                       0.161,
                       0.075,
                       0.000005)
calib.params [3,] <- c(paste0(data.path, cal.file.3),
                       84,
                       133,
                       0.087,
                       0.059,
                       0.000005)

### set up matrix to hold calibration results ###

calib.results <- matrix(ncol = 11, nrow = calib.n)
row.names(calib.results) <- c(1:calib.n)
colnames(calib.results) <- c('L0', 'T0', 'R2', 'k', 'half.life', 'start.buffer', 'end.buffer', 'start.std', 'end.std', 'start.sod', 'end.sod')

### calculate standard curve ###

c <- 2 ## testing

for(c in 1:calib.n){

  ### manually entered parameters ###
  
  spec.time <- as.numeric(calib.params[c, 2]) # time standard was read on spec
  sample.time <- as.numeric(calib.params[c, 3]) # time ko2 goes on the instrument
  calib.240 <- as.numeric(calib.params[c, 4]) - as.numeric(calib.params[c, 5]) # abs at 240
  aliquot <- as.numeric(calib.params[c, 6]) # size of standard aliquot in L
  
  ## read in times
  
  calib <- read.table(calib.params[c, 1], skip = 1)
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

## write so that you don't need to redo the calibration

write.table(calib.results, paste0(data.path, calib.file), quote = F, row.names = F)

#### calculate calibration factor ####
## if you don't want to redo the calibration
## you can start here

calib.results <- read.table(paste0(data.path, calib.file), header = T)

## if you need to drop a calibration point do it here
calib.results <- calib.results[-1,]

plot(c(calib.results[,2], 0) ~ c(calib.results[,1], 0),
     xlab = 'L0',
     ylab = '[SO] at t0 (pM)')

calib.lm <- lm(calib.results[,2] ~ 0 + calib.results[,1])
abline(calib.lm)
summary(calib.lm)

calib.factor <- as.numeric(calib.lm$coefficients[1])

#### assays ####

colnames(data.files) <- c('files')

## create matrix to hold output

so.matrix <- matrix(nrow = length(data.files$files), ncol = 36)
rownames(so.matrix) <- data.files$files
colnames(so.matrix) <- c(paste0(sample1, '.steady.state.conc'),
                         paste0(sample1, '.production.rate'),
                         paste0(sample1, '.cor.production.rate'),
                         paste0(sample2, '.steady.state.conc'),
                         paste0(sample2, '.production.rate'),
                         paste0(sample2, '.cor.production.rate'),
                         paste0(sample3, '.steady.state.conc'),
                         paste0(sample3, '.production.rate'),
                         paste0(sample3, '.cor.production.rate'),
                         paste0(sample4, '.steady.state.conc'),
                         paste0(sample4, '.production.rate'),
                         paste0(sample4, '.cor.production.rate'),
                         'buffer.start',
                         'buffer.end',
                         paste0(sample1, '.start'),
                         paste0(sample1, '.end'),
                         paste0(sample2, '.start'),
                         paste0(sample2, '.end'),
                         paste0(sample3, '.start'),
                         paste0(sample3, '.end'),
                         paste0(sample4, '.start'),
                         paste0(sample4, '.end'),
                         'sod.start',
                         'sod.end',
                         'buffer.mean',
                         'buffer.sd',
                         paste0(sample1, '.mean'),
                         paste0(sample1, '.sd'),
                         paste0(sample2, '.mean'),
                         paste0(sample2, '.sd'),
                         paste0(sample3, '.mean'),
                         paste0(sample3, '.sd'),
                         paste0(sample4, '.mean'),
                         paste0(sample4, '.sd'),
                         'sod.mean',
                         'sod.sd')

file <- data.files$files[3]

calc.so <- function(file){
  
  assay <- read.table(paste0(data.path, file), skip = 1)
  assay$V1 <- assay$V1 * 2 # convert from events to seconds
  colnames(assay) <- c('s', 'lum')
  
  plot(assay$lum ~ assay$s,
       type = 'l',
       ylab = 'Luminescence',
       xlab = 'Time (s)')
  
  ## read in times and plot
  
  print(paste(file, ':identify two stable buffer points'))
  buffer.start <- assay$s[identify(assay$s, assay$lum, n = 1, plot = T)]
  buffer.end <- assay$s[identify(assay$s, assay$lum, n = 1, plot = T)]
  
  print(paste(file, ':identify cell signal 1 start and end points'))
  sample1.start <- assay$s[identify(assay$s, assay$lum, n = 1, plot = T)]
  sample1.end <- assay$s[identify(assay$s, assay$lum, n = 1, plot = T)]
  
  print(paste(file, ':identify cell signal 2 start and end points'))
  sample2.start <- assay$s[identify(assay$s, assay$lum, n = 1, plot = T)]
  sample2.end <- assay$s[identify(assay$s, assay$lum, n = 1, plot = T)]
  
  print(paste(file, ':identify cell signal 3 start and end points'))
  sample3.start <- assay$s[identify(assay$s, assay$lum, n = 1, plot = T)]
  sample3.end <- assay$s[identify(assay$s, assay$lum, n = 1, plot = T)]
  
  print(paste(file, ':identify cell signal 4 start and end points'))
  sample4.start <- assay$s[identify(assay$s, assay$lum, n = 1, plot = T)]
  sample4.end <- assay$s[identify(assay$s, assay$lum, n = 1, plot = T)]
  
  print(paste(file, 'identify two points along post-SOD baseline'))
  sod.start <- assay$s[identify(assay$s, assay$lum, n = 1, plot = T)]
  sod.end <- assay$s[identify(assay$s, assay$lum, n = 1, plot = T)] 
  
  ## average values for times
  
  buffer.mean <- mean(assay$lum[assay$s >= buffer.start & assay$s <= buffer.end])
  buffer.sd <- sd(assay$lum[assay$s >= buffer.start & assay$s <= buffer.end])
  
  sample1.mean <- mean(assay$lum[assay$s >= sample1.start & assay$s <= sample1.end])
  sample1.sd <- sd(assay$lum[assay$s >= sample1.start & assay$s <= sample1.end])
  
  sample2.mean <- mean(assay$lum[assay$s >= sample2.start & assay$s <= sample2.end])
  sample2.sd <- sd(assay$lum[assay$s >= sample2.start & assay$s <= sample2.end])
  
  sample3.mean <- mean(assay$lum[assay$s >= sample3.start & assay$s <= sample3.end])
  sample3.sd <- sd(assay$lum[assay$s >= sample3.start & assay$s <= sample3.end])
  
  sample4.mean <- mean(assay$lum[assay$s >= sample4.start & assay$s <= sample4.end])
  sample4.sd <- sd(assay$lum[assay$s >= sample4.start & assay$s <= sample4.end])
  
  sod.mean <- mean(assay$lum[assay$s >= sod.start & assay$s <= sod.end])
  sod.sd <- sd(assay$lum[assay$s >= sod.start & assay$s <= sod.end])
  
  ## calculate final values
  
  sample.sig1 <- sample1.mean - buffer.mean
  
  steady.state.so1 <- (sample.sig1 - (buffer.mean - sod.mean)) * calib.factor # steady state so concentration in pM
  net.prod.so1 <- steady.state.so1 * flow.rate / 1000 # net so production in pmol ml-1 min-1
  cor.prod.so1 <- net.prod.so1 / calib.factor # so production corrected for decay in pmol ml-1 min-1
  
  sample.sig2 <- sample2.mean - buffer.mean
  
  steady.state.so2 <- (sample.sig2 - (buffer.mean - sod.mean)) * calib.factor # steady state so concentration in pM
  net.prod.so2 <- steady.state.so2 * flow.rate / 1000 # net so production in pmol ml-1 min-1
  cor.prod.so2 <- net.prod.so2 / calib.factor # so production corrected for decay in pmol ml-1 min-1
  
  sample.sig3 <- sample3.mean - buffer.mean
  
  steady.state.so3 <- (sample.sig3 - (buffer.mean - sod.mean)) * calib.factor # steady state so concentration in pM
  net.prod.so3 <- steady.state.so3 * flow.rate / 1000 # net so production in pmol ml-1 min-1
  cor.prod.so3 <- net.prod.so3 / calib.factor # so production corrected for decay in pmol ml-1 min-1
  
  sample.sig4 <- sample4.mean - buffer.mean
  
  steady.state.so4 <- (sample.sig4 - (buffer.mean - sod.mean)) * calib.factor # steady state so concentration in pM
  net.prod.so4 <- steady.state.so4 * flow.rate / 1000 # net so production in pmol ml-1 min-1
  cor.prod.so4 <- net.prod.so4 / calib.factor # so production corrected for decay in pmol ml-1 min-1
  
  ## generate a plot illustrating calculations
  
  pdf(paste0(data.path, file, '.pdf'), width = 6, height = 5)
  
  plot(assay$lum ~ assay$s,
       type = 'l',
       ylab = 'Luminescence',
       xlab = 'Time (s)')
  
  lines(c(buffer.start, buffer.start),
        c(-100, 1e6),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(buffer.end, buffer.end),
        c(-100, 1e6),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sample1.start, sample1.start),
        c(-100, 1e6),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sample1.end, sample1.end),
        c(-100, 1e6),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sample2.start, sample2.start),
        c(-100, 1e6),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sample2.end, sample2.end),
        c(-100, 1e6),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sample3.start, sample3.start),
        c(-100, 1e6),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sample3.end, sample3.end),
        c(-100, 1e6),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sample4.start, sample4.start),
        c(-100, 1e6),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sample4.end, sample4.end),
        c(-100, 1e6),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sod.start, sod.start),
        c(-100, 1e6),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sod.end, sod.end),
        c(-100, 1e6),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(buffer.start, buffer.end),
        c(buffer.mean, buffer.mean),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sample1.start, sample1.end),
        c(sample1.mean, sample1.mean),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sample2.start, sample2.end),
        c(sample2.mean, sample2.mean),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sample3.start, sample3.end),
        c(sample3.mean, sample3.mean),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sample4.start, sample4.end),
        c(sample4.mean, sample4.mean),
        lwd = 0.5,
        col = 'grey')
  
  lines(c(sod.start, sod.end),
        c(sod.mean, sod.mean),
        lwd = 0.5,
        col = 'grey')
  
  dev.off()
   
  print(c(steady.state.so1, net.prod.so1, cor.prod.so1))
  print(c(steady.state.so2, net.prod.so2, cor.prod.so2))
  print(c(steady.state.so3, net.prod.so3, cor.prod.so3))
  print(c(steady.state.so4, net.prod.so4, cor.prod.so4))
  
  return(c(steady.state.so1,
           net.prod.so1,
           cor.prod.so1,
           steady.state.so2,
           net.prod.so2,
           cor.prod.so2,
           steady.state.so3,
           net.prod.so3,
           cor.prod.so3,
           steady.state.so4,
           net.prod.so4,
           cor.prod.so4,
           buffer.start,
           buffer.end,
           sample1.start,
           sample1.end,
           sample2.start,
           sample2.end,
           sample3.start,
           sample3.end,
           sample4.start,
           sample4.end,
           sod.start,
           sod.end,
           buffer.mean,
           buffer.sd,
           sample1.mean,
           sample1.sd,
           sample2.mean,
           sample2.sd,
           sample3.mean,
           sample3.sd,
           sample4.mean,
           sample4.sd,
           sod.mean,
           sod.sd))  
}

for(file in data.files$files){  
  so.matrix[which(row.names(so.matrix) == file),] <- calc.so(file)
}

write.table(so.matrix, paste0(data.path, assay.file), quote=F)
