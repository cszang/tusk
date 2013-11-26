## tusk is small library to help with commonly used gridded data sets.
## ATM, the following data sets are supported:
## - CRU TS 3.20

##' Reported supported data sets
##'
##' ts320: CRU TS 3.20
##' @title Supported data sets
##' @return a vector holding the short names of supported data sets
##' @examples
##' supported_sets()
##' @export
supported_sets <- function() {
  supp <- list(abbrev = c("ts320"),
       full = c("CRU TS 3.20"))
  class(supp) <- c("list", "tusk_supported_sets")
  supp
}

##' @export
print.tusk_supported_sets <- function(x, ...) {
  n <- length(x$abbrev)
  m <- matrix(NA, nrow = n, ncol = 2)
  for (i in 1:n) {
    m[i,] <- c(x$abbrev[i], x$full[i])
  }
  m <- data.frame(Abbreviation = m[,1], Full.Name = m[,2])
  print(m, ...)
}

##' Convert a coordinate from decimal degrees to radians
##' 
##' This function converts coordinates from decimal degrees to radians.
##' @title Convert degrees to radians
##' @param x the coordinate in decimal degrees
##' @return the coordinate as radians
##' @examples
##' degrees_to_radians(40.81)
##' @export
degrees_to_radians <- function(x) {
  degrees <- (x %/% 1)
  minutes <- x - degrees
  decimal <- minutes/60
  c.num <- degrees + decimal
  radians <- c.num * pi / 180
  radians
}

##' Calculate the distance between two points as great circle distance
##' on the globe
##'
##' Uses information from
##' http://janmatuschek.de/LatitudeLongitudeBoundingCoordinates, where
##' dist = arccos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) *
##' cos(lon1 - lon2)) * R
##' @title Great circle distance
##' @param x1 coordinate pair 1, as list with $lat and $lon
##' @param x2 coordinate pair 2, as list with $lat and $lon
##' @return The distance on the globe in kilometers.
##' @examples
##' great_circle_dist(list(lon = 0.3, lat = 40.91), list(lon = 1.3,
##' lat = 41.02))
##' @export
great_circle_dist <- function(x1, x2) {
  X1 <- degrees_to_radians(x1$lat)
  X2 <- degrees_to_radians(x2$lat)
  Y1 <- degrees_to_radians(x1$lon)
  Y2 <- degrees_to_radians(x2$lon)
  dist <- acos(sin(X1) * sin(X2) + cos(X1) * cos(X2) * cos(Y1 - Y2)) * 6370
  dist
}

##' Get four nearest gridpoints for a lon/lat pair.
##'
##' The returned list of gridpoint coordinates can be used for
##' interpolation
##' @title Get four nearest gridpoints
##' @param coords the coordinates to find the nearest gridpoints for
##' as list with $lon and $lat in decimal degrees
##' @param netcdf the ncdf object holding the relevant data
##' @param data_set the kind of data used (see ?supported_sets) for details
##' @return a list holding four lists, with each holding a lat/lot
##' @import ncdf
##' @export
four_nearest <- function(coords, netcdf, data_set = "ts320") {
  if (!any(supported_sets()$abbrev == data_set)) {
    stop("Data set not supported. See ?supported_sets.")
  }
  require(ncdf)
  if (data_set == "ts320") {
    all_lon <- get.var.ncdf(netcdf, "lon")
    all_lat <- get.var.ncdf(netcdf, "lat")
  }
  lon <- coords$lon
  lat <- coords$lat
  lons <- c(tail(which(all_lon < lon), 1), which(all_lon > lon)[1])
  lats <- c(tail(which(all_lat < lat), 1), which(all_lat > lat)[1])
  # combine to four gridpoints
  list(
    list(lon = lons[1], lat = lats[1]),
    list(lon = lons[1], lat = lats[2]),
    list(lon = lons[2], lat = lats[2]),
    list(lon = lons[2], lat = lats[1])
  )
}

##' Interpolate gridded climate data for a given coordinate pair
##'
##' Interpolation is done using inverse distance weighting over the
##' four nearest grid points. The results will depend on the kind of
##' data used for interpolation. E.g., using CRU TS 3.20 will result
##' in monthly data starting in January 1901.
##' @title Interpolate from gridded climate data sets
##' @param netcdf the ncdf object the data should be extracted from
##' @param param the name of the paramater to extract as character
##' string
##' @param coords the coordinates of the point to get the data for as
##' list with $lon and $lat
##' @param nearest a nearest gridpoints objects as returned from
##' four_nearest()
##' @param data_set the kind of data set used (see ?supported_sets)
##' @return a data.frame holding the extracted and interpolated data
##' @examples
##' \dontrun{
##' library(ncdf)
##' cru_maxtemp <- open.ncdf("~/Data/cru_ts3.20.1901.2011.tmx.dat.nc")
##' my_coords <- list(lon = 0.3, lat = 40.81)
##' my_nearest <- four_nearest(my_coords, cru_maxtemp, "ts320")
##' interpol_four(cru_maxtemp, "tmx", my_coords, my_nearest, "ts320")
##' }
##' @import ncdf
##' @export
interpol_four <- function(netcdf, param, coords, nearest, data_set = "ts320") {
  if (!any(supported_sets()$abbrev == data_set)) {
    stop("Data set not supported. See ?supported_sets.")
  }
  require(ncdf)
  dists <- numeric(4)
  extract <- NULL
  for (i in 1:4) {
    this_lon <- nearest[[i]]$lon
    this_lat <- nearest[[i]]$lat
    .start <- switch(data_set,
                     ts320 = c(this_lon, this_lat, 1))
    .count <- switch(data_set,
                     ts320 = c(1, 1, -1))
    extract <- cbind(extract,
                     get.var.ncdf(netcdf, param,
                                  start = .start,
                                  count = .count))
    dists[i] <- great_circle_dist(nearest[[i]],
                                  list(lon = coords$lon,
                                       lat = coords$lat))
  }
  weights <- 1/dists
  extract_weight <- sweep(extract, 2, weights, `*`)
  extract_int <- apply(extract_weight, 1, function(x) sum(x/sum(weights)))
  
  if (data_set == "ts320") {
    first_date <- as.Date("1901-01-01")
    whole_period <- seq(first_date, length.out = 1332, by = "1 month")
    Years <- as.numeric(substr(whole_period, 1, 4))
    Months <- as.numeric(substr(whole_period, 6, 7))
    out <- data.frame(
      year = Years,
      month = Months,
      extract = extract_int
    )
    out
  }
}

