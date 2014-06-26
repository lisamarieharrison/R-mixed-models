#input date time as a single character
a <- "20060121 073819"

#change time to POSIX* class
b <- strptime(a, format = "%Y%m%d %H%M%S")

hrs <- function(u) {
  x <- u * 3600
  return(x)
}

mns <- function(m) {
  x <- m * 60
  return(x)
}


