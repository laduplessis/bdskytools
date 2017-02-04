
dateorigin <- "1970-01-01"

#' Convert date from year + fraction to normal date (day, month year)
#' 
#' Does not work for years with no decimal (should be December 31 of the previous year)
#' 
#' @export
getDayDate <- function(date) {

  year <- floor(date)
  len  <- as.numeric(format.Date(as.Date(paste0(year,"-12-31"),),"%j"))
  day  <- floor((date-year)*len)
  
  #return(as.Date(paste(year,day), format="%Y %j"))
  return(as.Date(paste0(year-1,"-12-31"))+day)
}

#' Convert from date to decimal fraction of year
#' 
#' @export
getYearDate <- function(date) {
  
  d    <- as.Date(date)
  year <- as.numeric(format.Date(d,"%Y"))
  day  <- as.numeric(format.Date(d,"%j"))
  len  <- as.numeric(format.Date(as.Date(paste0(year,"-12-31"),),"%j"))
  
  return (year + day/len)
}

#' From the first week BEFORE start up to end (strict)
#' 
#' Weeks start on Sunday
#' 
#' @param strictend: If TRUE the last value is end, if FALSE it is the first week after end.
#' @param inclusive: If TRUE then the weeks start and end m
#' 
#' @export
getWeeks <- function(start, end, inclusive=TRUE, strictend=FALSE) {
  
  weeks <- seq(from=as.numeric(start)-as.numeric(format.Date(start, format="%u")), to=as.numeric(end), by=7)
  if (max(weeks) < end) {
    if (strictend == TRUE)
      weeks <- c(weeks, end)
    else { 
      if (inclusive == TRUE)
          weeks <- c(weeks, max(weeks)+7)
    }
  }
  
  if (inclusive == FALSE) {
    weeks <- weeks[2:length(weeks)]
  }
  
  return(as.Date(weeks, origin=dateorigin))
}


#' Return start and end date, with all months in between
#' 
#' From start to end, breaking at the 1st of every month in-between
#' 
#' @export
getMonths <- function(start, end) {
  
  startmonth <- as.numeric(format.Date(start,"%m"))+1
  endmonth   <- as.numeric(format.Date(end,"%m"))
  startyear  <- as.numeric(format.Date(start,"%Y"))
  endyear    <- as.numeric(format.Date(end,"%Y"))
  
  months <- c(start)
  for (year in startyear:endyear) {
    
    if (year == startyear) 
      yearstart = startmonth
    else
      yearstart = 1
    
    if (year == endyear)
      yearend = endmonth
    else
      yearend = 12
    
    if (yearstart <= yearend) {
      for (m in yearstart:yearend) {
        months <- c(months, as.Date(sprintf("%d-%d-1",year,m)))
      }
    }    
  }
  months <- c(months, end)
  
  return(months)
  
}