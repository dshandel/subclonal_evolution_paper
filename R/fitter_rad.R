fitter <- function(df) {
  
  # defining variables to work with
  d <- df
  
  
  lines <- character()
  # for loop, checks if CF$cline is in vector lines, if so, do nothing, else append to vector. 
  for (line in d$cline) (
    ifelse(line %in% lines, NA, lines <- c(lines, line)))
  
  # separate into dataframes according to cline
  per_line <- split(d, f = d$cline, drop = TRUE)
  
  
  # Function to make model for every line and dataframe with predicted data for each line. 
  fit_model= function(line) {
    # line is a dataframe. Outputs the fits of norm and dose
    fit <- drm(relative ~ dose, # define y -axis (ncolonies) and x-axis (dose)
               data = line, # defines dataframe
               fct = LL.4 (names = c('Slope', "Lower Limit", "Upper Limit", "IC50"))) # defines to fit a log-losistic model)
    
    newdata <- expand.grid(dose=exp(seq(log(0.00000001), log(max(df$dose) + 1000), length=1000))) # new data with doses. Note: lowest dose is not
    # log 0 but log('very small number') because otherwise this will hamper the scaling in ggplot later on. 
    pm <- predict(fit, newdata=newdata, interval="confidence") # new data with predictions and confidence intervals
    newdata$pred <- pm[,1] # add prediction values to new data.
    newdata$predmin <- pm[,2] # add lower bounderies to new data 
    newdata$predmax <- pm[,3] # add upper bounderies to new data 
    newdata$cline <- line$cline[1] # add column with cline.
    newdata$expcode = unique(df$expcode)
    return(newdata)
  }
  
  # Store in vector
  data_frames <- lapply(per_line, fit_model) 
  
  
  #merge dataframes
  df = do.call(rbind, data_frames)
  
  return(df) 
  
} 