library(dplyr)
library(seewave) # for bwfilter https://www.rdocumentation.org/packages/seewave/versions/2.0.5/topics/bwfilter 
library(pracma) # for findpeaks https://www.rdocumentation.org/packages/pracma/versions/1.9.9/topics/findpeaks

peak_data <- read.csv("peakdata.csv")

# filter settings
padlength <- 100 # number of samples to pad for better filtering
sampling_freq <- 1000
cutoff <- 15 # remove frequencies above this

# settings for peak picking
nups <- 10 # minimum number of increasing steps before a peak is reached
ndowns <- 10 # minimum number of decreasing steps after the peak
minpeakheight <- 1.5 # the minimum (absolute) height a peak has to have to be recognized as such
threshold <- 2 # the minimum value to be classified as peak
max_wave_length = 250 # maximum number of samples to seek for a trough after peak
largest = FALSE # true for highest peak, false for first peak

# settings for trough
minpeak <- 0.5 # tough must be this much lower than the peak to be classified


p2p_amp <- function(time, data, visualise=FALSE, filename=""){
  
  # reorder in time
  order_index <- order(time)
  data <- data[order_index]
  
  # first pad the end of the data to stop filter dipping at end
  len <- length(data)
  data[len : (len + padlength)] <- data[len]
  # smooth data
  filt_data <- seewave::bwfilter(data, f=sampling_freq, n=2, to=cutoff) %>% as.vector

  # delete pad we added
  data <- data[1:len]
  filt_data <- filt_data[1:len]

  # get first peak meeting criteria
  peak <- pracma::findpeaks(filt_data,
                            npeaks = 1,
                            nups = nups,
                            ndowns = ndowns,
                            minpeakheight = minpeakheight,
                            sortstr=largest)
  
  # if we found a suitable peak, get value of peak, otherwise stop and return NA
  if (!is.null(peak)){
    max_val = peak[1]
    max_val_index = peak[2]
  } else {
    if (visualise) plot_orig_and_smooth(time, data, filt_data, "No peak")
    warning(sprintf("! no peak found for %s\n", filename))
    return(NA)
  }
  
  data_after_peak <- filt_data[max_val_index : (max_val_index + max_wave_length)] # data after peak up to max wave length or end of data
  data_after_peak <- data_after_peak[!is.na(data_after_peak)]
  
  
  # get first trough meeting criteria
  min_trough_val <- max_val - minpeak
  trough <- pracma::findpeaks(-1 * data_after_peak, # negate to find troughs
                              npeaks = 1,
                              nups = nups,
                              ndowns = ndowns,
                              minpeakheight = -1 * min_trough_val,
                              sortstr=largest)
  
  # if we found a suitable trough, get value of trough, otherwise stop and return NA
  # we could add a fallback to just get minimum value rather than a trough
  if (!is.null(trough)){
    min_val = -1 * trough[1]
    min_val_index = max_val_index + trough[2] - 1
  } else {
    if (visualise){
      plot_orig_and_smooth(time, data, filt_data, "No trough")
      points(c(time[max_val_index]), c(max_val), pch=0, col=30, cex=3) # square
    }
    warning(sprintf("! no trough found for %s\n", filename))
    return(NA)
  }

  # peak to peak amplitude
  p2p <- max_val - min_val
  
  # visualisation - plot graphs to inspect
  if (visualise){
    
    plot_orig_and_smooth(time, data, filt_data, "OK")
    points(c(time[max_val_index], time[min_val_index]), c(max_val, min_val), pch=0, col=30, cex=3) # squares
    
  }
  
  return(p2p)
}

# plot function
plot_orig_and_smooth <- function(time, x, xsmooth, ttl=""){
  plot(time, x, xlim=c(100,400), ylim=c(-4,12)) # points
  lines(time, xsmooth) # smoothed line
  title(main=ttl)
}


#############################################

amplitudes <- peak_data %>%
  dplyr::group_by(File,condition) %>%
  filter(variable =="FCZ") %>% 
  summarise(p2p_amp = p2p_amp(Time, value, visualise = TRUE, filename = File[1])) # run function on groups

