IC50_fun <- function(df, rad) {
  
  # defining variables to work with
  ifelse(is.null(rad), d<-df, d <- subset(df, rec_rad == rad))
  
  lines <- character()
  # for loop, checks if CF$cline is in vector lines, if so, do nothing, else append to vector. 
  for (line in d$cline) (
    ifelse(line %in% lines, NA, lines <- c(lines, line)))
  
  # separate into dataframes according to cline
  per_line <- split(d, f = d$cline, drop =TRUE)
  
  
  # Function to make model for every line and dataframe with predicted data for each line. 
  ic50= function(line) {
    # line is a dataframe. Outputs the fits of norm and dose
    fit <- drm(relative ~ dose, # define y -axis (ncolonies) and x-axis (dose)
               data = line, # defines dataframe
               fct = LL.4 (names = c('Slope', "Lower Limit", "Upper Limit", "IC50"))) # defines to fit a log-losistic model)
    
    return(coef(fit)[4])
  }
  
  # Store in vector
  IC50_vector <- lapply(per_line, ic50) 
  
  #merge dataframes
  df = do.call(rbind, IC50_vector)
  
  return(df) 
  
} 