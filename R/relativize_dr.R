relative <- function(df, rad) {
  #' @param df with dose response data and columns cline, dose, Exp, 
  #' rec_rad, expcode, ncolonies
  #' @param rad if line received radiation (T or F), stored in column rec_rad
  #' @return relativized data where dose = 0 is used as 100%
  
  # defining variables to work with
  ifelse(rad == T, d<-subset(df, rec_rad == 1), d<-subset(df, rec_rad == 0) )
  
  # define ncolonies
  column <- "ncolonies"
  
  # Extracts the column name.
  col <- d[c(column)]
  
  # extracts number of experiemnts 
  n <- length(unique(d$Exp))
  
  # Takes mean every nth row
  
  mean_fun <- function(x) {
    m <- mean(x, na.rm = TRUE)
    return(m)
  }
  
  mean <- aggregate(col,list(rep(1:(nrow(col)%/%n+1),each=n,len=nrow(col))), mean_fun)[-1]
  # Duplicates (n=2) or triplicates (n=3) the rows.
  mean <- mean[rep(seq_len(nrow(mean)), each = n), ]
  # Adds everything to the dataframe. 
  d["mean"] <- mean
  
  # Select mean values of 0 concentration
  d <- d %>%
    group_by(cline) %>%
    arrange(cline,dose)
  
  first <- d[d$dose==0, ] %>%
    dplyr::select(cline, value100=mean)
  
  # Make cline a factor
  d$cline<-factor(d$cline)
  
  # Merge first original (d)
  d <- d %>%
    merge(first,by=c("cline"))
  
  # Extract only every nth row.
  d = d[seq(1, nrow(d), n), ]
  
  # make new column with relative, called relative
  d$relative = d$ncolonies/d$value100
  d$relative_mean = d$mean/d$value100
  
  d<-subset(d, select=-c(mean,value100))
  
  # compute relative standard error of the mean
  # define the column with sem in it
  column <- "relative"
  
  # Extracts the column name.
  col <- d[c(column)]
  
  
  sem_fun <- function(x) {
    std <- sd(x, na.rm = TRUE)
    vector <- na.omit(x)
    
    sem <- std/sqrt(length(vector))
    return(sem)
  }
  
  sem <- aggregate(col,list(rep(1:(nrow(col)%/%n+1),each=n,len=nrow(col))),sem_fun)[-1]
  sem <- sem[rep(seq_len(nrow(sem)), each = n), ]
  d["relative_sem"] <- sem
  
  return(d)
  
}