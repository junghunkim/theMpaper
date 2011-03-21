w18w38w58.map <- function(str){
  w18.start.time <- as.POSIXct("2000-12-18 01:00:00 GMT", tz = "GMT")
  w38.start.time <- as.POSIXct("2001-05-07 01:00:00 GMT", tz = "GMT")
  w58.start.time <- as.POSIXct("2001-09-24 01:00:00 GMT", tz = "GMT")
  w78.start.time <- as.POSIXct("2002-02-11 01:00:00 GMT", tz = "GMT")

  d <- as.POSIXct(as.character(str), "%a, %d %b %Y %H", tz = "GMT")

  if(w18.start.time <= d && d < w38.start.time)
    return("w18")
  else if(w38.start.time <= d && d < w58.start.time)
    return("w38")
  else if(w58.start.time <= d && d < w78.start.time)
    return("w58")
  else
    return("NA")
}

unit.interval <- function(date, labels.vec, label){

  idx <- which(labels.vec == label)
  unit.interval.vec <- seq(0,0,length.out = length(idx))
  
  date.min <- as.POSIXct(as.character(date[idx[1]]), "%a, %d %b %Y %H", tz = "GMT")
  date.max <- as.POSIXct(as.character(date[idx[length(idx)]]),
                         "%a, %d %b %Y %H", tz = "GMT")

  duration <- as.numeric(date.max - date.min, unit = "hours")
  
  for(i in 1:length(idx)){
    
    date.i <- as.POSIXct(as.character(date[idx[i]]), "%a, %d %b %Y %H", tz = "GMT")
    unit.interval.vec[i] = as.numeric((date.i - date.min), unit = "hours")/duration
  }
  return(unit.interval.vec)
}
