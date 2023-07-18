sci_not <- function(n, round_n) {
  #' Converts a number to scientific notation and extracts number and power
  #' @param n a number to convert
  #' @param round_n how to round the number
  #' @return a list with the number and the power
  require(stringr)
  
  output <- format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output <- sub("e", "*10^", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  
  number <- round(as.numeric(str_split(output, pattern = '\\*')[[1]][1]), round_n)
  power <- (str_split(output, pattern = '\\^')[[1]][2])
  
  return(c(number, power))
  
}


n = mwu$p
round_n = 2
