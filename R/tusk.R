##' Reported supported data sets
##'
##' ts322: CRU TS 3.22
##' ts321: CRU TS 3.21
##' spei22: SPEIbase v.2.2
##' pdsidai2011: scPDSI from Dai 2011a and 2011b
##' eobs140: E-OBS data daily data version 14.0
##' puhg_pet: Princeton University Hydroclimatology Group 61-yr (1948-2008) Potential Evaporation Dataset
##' @title Supported data sets
##' @return a vector holding the short names of supported data sets
##' @examples
##' supported_sets()
##' @export
supported_sets <- function() {
  supp <- list(abbrev = c("ts322", "ts321", "spei22", "pdsidai2011",
                         "eobs140", "puhg_pet"),
       full = c("CRU TS 3.22", "CRU TS 3.21", "SPEIbase v.2.2",
                "scPDSI from Dai 2011a and 2011b", "E-OBS 14.0",
                "Princeton University Hydroclimatology Group 61-yr (1948-2008) Potential Evaporation Dataset"))
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
  x * pi / 180
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

##' Get four or sixteen nearest gridpoints for a lon/lat pair.
##'
##' The returned list of gridpoint coordinates can be used for
##' interpolation
##' @title Get 4 or 16 nearest gridpoints
##' @param coords the coordinates to find the nearest gridpoints for
##' as list with $lon and $lat in decimal degrees
##' @param netcdf the ncdf object holding the relevant data
##' @param data_set the kind of data used (see ?supported_sets) for details
##' @return a list holding four or sixteen lists, with each holding a
##' lat/lot
##' @import ncdf4
##' @export
nearest_points <- function(coords, netcdf, data_set = "ts322", npoints = 4) {
  if (!any(supported_sets()$abbrev == data_set)) {
    stop("Data set not supported. See ?supported_sets.")
  }
  if (any(data_set == c("ts322", "ts321", "spei22", "pdsidai2011"))) {
    all_lon <- ncvar_get(netcdf, "lon")
    all_lat <- ncvar_get(netcdf, "lat")
  } else {
    if (any(data_set == c("eobs140", "puhg_pet"))) {
      all_lon <- ncvar_get(netcdf, "longitude")
      all_lat <- ncvar_get(netcdf, "latitude")
    }
  }
  lon <- coords$lon
  lat <- coords$lat
  lons <- c(tail(which(all_lon < lon), 1), which(all_lon > lon)[1])
  lats <- c(tail(which(all_lat < lat), 1), which(all_lat > lat)[1])
  replace_missing <- function(x) {
    if (length(x) == 1) {
      x <- c(x, x)
    }
    x
  }
  lons <- replace_missing(lons)
  lats <- replace_missing(lats)
  lons16 <- c(min(lons) - 1, min(lons), max(lons), max(lons) + 1)
  lats16 <- c(min(lats) - 1, min(lats), max(lats), max(lats) + 1)

  if (npoints == 4) {
    gridm <- expand.grid(lon = lons, lat = lats)
  } else {
    gridm <- expand.grid(lon = lons16, lat = lats16)
  }

  index <- degrees <- list()
  o <- dim(gridm)[1]

  for (i in 1:o) {
    index[[i]] <- list(lon = gridm$lon[i], lat = gridm$lat[i])
    degrees[[i]] <- list(lon = all_lon[gridm$lon[i]],
                         lat = all_lat[gridm$lat[i]])
  }

  class(index) <- c("tusk_nearest_points", "list")
  attributes(index)$coords <- degrees
  index
}

##' @export
print.tusk_nearest_points <- function(x, ...) {
  n <- length(x)
  coords <- attributes(x)$coords
  pretty_coords <- matrix(NA, ncol = 2, nrow = n)
    for (i in 1:n) {
    pretty_coords[i,] <- c(coords[[i]]$lon, coords[[i]]$lat)
  }
  pretty_coords <- data.frame(pretty_coords)
  names(pretty_coords) <- c("Longitude", "Latitude")
  print(pretty_coords, ...)
}

##' Interpolate and downscale gridded climate data for a given
##' coordinate pair
##'
##' Interpolation is done using inverse distance weighting over the
##' four nearest grid points. The results will depend on the kind of
##' data used for interpolation. E.g., using CRU TS 3.22 will result
##' in monthly data starting in January 1901. For CRU precipitation
##' and temperature data, the data can be downscaled to 30 arc seconds
##' resolution using the Worldclim climatology. In this case, for the
##' given points, the anomalies of the nearest gridpoints are
##' calculated based on their local 1950-2000 climatology, and
##' interpolated. The interpolated anomalies are then rescaled to the
##' Worldclim climatology for the nearest point in the Worldclim grid.
##' @title Interpolate and downscale from gridded climate data sets
##' @param netcdf the ncdf object the data should be extracted from
##' @param worldclim path to folder with worldclim .bil files
##' (optional, and only used if \code{downscale} is TRUE)
##' @param param the name of the paramater to extract as character
##' string
##' @param coords the coordinates of the point to get the data for as
##' list with $lon and $lat
##' @param nearest a nearest gridpoints objects as returned from
##' four_nearest()
##' @param data_set the kind of data set used (see ?supported_sets)
##' @param downscale logical: shall the data be downscaled to 30 arc
##' seconds?
##' @return a data.frame holding the extracted and interpolated data
##' @examples
##' \dontrun{
##' library(ncdf4)
##' cru_maxtemp <- nc_open("~/Data/cru_ts3.22.1901.2013.tmp.dat.nc")
##' my_coords <- list(lon = 0.3, lat = 40.81)
##' my_nearest <- nearest_points(my_coords, cru_maxtemp, "ts322")
##' interp_down(cru_maxtemp, "tmx", my_coords, my_nearest, "ts322")
##' }
##' @import ncdf4 raster rgdal
##' @export
interp_down <- function(netcdf, worldclim = NULL, param, coords,
                             nearest, data_set = "ts322", downscale = FALSE) {

  if (!any(supported_sets()$abbrev == data_set)) {
    stop("Data set not supported. See ?supported_sets.")
  }

  if (downscale & is.null(worldclim)) {
    stop("For downscaling worldclim data is needed.")
  }

  if (downscale & (!any(c("ts322", "ts321") == data_set))) {
    stop("Downscaling only implemented for CRU data.")
  }

  ## methods for computing anomalies (meth1) and rescaling to
  ## worldclim climatology (meth2): this is subtraction/addition for
  ## temperatures, and division/multiplication for precipitation
  if (downscale) {
    if(any(c("prec", "pre") == param)) {
      meth1 <- function(x, y) { y / x }
      meth2 <- function(x, y) { y * x }
      scale_fac <- 1 # no division needed, WC is "as is"
    } else {
      meth1 <- function(x, y) { y - x }
      meth2 <- function(x, y) { y + x }
      scale_fac <- 10 # WC is T [degrees C] * 10
    }
  }

  first_date <- switch(data_set,
                      ts322 = as.Date("1901-01-01"),
                      ts321 = as.Date("1901-01-01"),
                      spei22 = as.Date("1901-01-01"),
                      pdsidai2011 = as.Date("1850-01-01"),
                      eobs140 = as.Date("1950-01-01"),
                      puhg_pet = as.Date("1948-01-01"))

  ## all dates in the data
  time_length <- length(ncvar_get(netcdf, "time"))
  if (any(c("eobs140") == data_set)) {
    dates_all <- seq(first_date, length.out = time_length,
                     by = "1 day")
  } else {
    dates_all <- seq(first_date, length.out = time_length,
                     by = "1 month")
  }
  ## the dates for the climatology reference period (according to
  ## worldclim)
  dates_clim <- seq(as.Date("1950-01-01"), as.Date("2000-12-31"),
                    by = "1 month")

  ## here, the selection for the reference climatoloy in the data
  ## start and end
  clim_start <- which(dates_all == dates_clim[1])
  clim_end <- which(dates_all == tail(dates_clim, 1))

  month_num <- 1:12
  ## format month number strings to have trailing 0 for months 1:9
  month_num_string <- sapply(month_num, function(x) {
    X <- as.character(x)
    if (nchar(X) < 2) {
      X <- paste("0", X, sep = "")
    }
    X
  })

  npoints <- length(nearest)

  dists <- numeric(npoints)
  extract <- NULL

  for (i in 1:npoints) {

    this_lon <- nearest[[i]]$lon
    this_lat <- nearest[[i]]$lat

    .start <- switch(data_set,
                    ts322 = c(this_lon, this_lat, 1),
                    ts321 = c(this_lon, this_lat, 1),
                    spei22 = c(this_lon, this_lat, 1),
                    pdsidai2011 = c(this_lon, this_lat, 1),
                    eobs140 = c(this_lon, this_lat, 1),
                    puhg_pet = c(this_lon, this_lat, 1, 1))

    .count <- switch(data_set,
                    ts322 = c(1, 1, -1),
                    ts321 = c(1, 1, -1),
                    spei22 = c(1, 1, -1),
                    pdsidai2011 = c(1, 1, -1),
                    eobs140 = c(1, 1, -1),
                    puhg_pet = c(1, 1, 1, -1))
    
    .extract <- ncvar_get(netcdf, param,
                         start = .start,
                         count = .count)

    if (downscale) {

      .extract_climperiod <- .extract[clim_start:clim_end]

      ## build climatologies for all months, compute anomalies

      for (j in month_num) {
        mreg <- paste(".*\\-", month_num_string[j], "\\-.*", sep = "")
        mmatch_clim <- grep(mreg, as.character(dates_clim))
        mmatch_all <- grep(mreg, as.character(dates_all))
        mmatch_clim <- mmatch_all[dates_all[mmatch_all] %in% dates_clim[mmatch_clim]]
        month_cell_clim <- .extract[mmatch_clim]
        month_cell_clim <- mean(month_cell_clim, na.rm = TRUE)
        .extract[mmatch_all] <- meth1(month_cell_clim, .extract[mmatch_all])
      }
    }

    extract <- cbind(extract, .extract)

    dists[i] <- great_circle_dist(attr(nearest, "coords")[[i]],
                                  list(lon = coords$lon,
                                       lat = coords$lat))
  }

  ## any points all NA (sea?)
  ##:ess-bp-start::browser@nil:##
  seap <- apply(extract, 2, function(x) all(is.na(x)))

  ## interpolate values/anomalies using inverse distance weighting
  weights <- 1/dists[!seap]
  if (sum(!seap) > 1) {
    extract_weight <- sweep(extract[,!seap], 2, weights, `*`)
    extract_int <- apply(extract_weight, 1, function(x)
      sum(x/sum(weights)))
  } else {
    if (sum(!seap) == 1) {
      extract_int <- extract[,!seap]
    } else {
      stop("All data is NA for the nearest points.")
    }
  }

  Years <- as.numeric(substr(dates_all, 1, 4))
  Months <- as.numeric(substr(dates_all, 6, 7))
  
  if (any(c("eobs140") == data_set)) {
    Days <- as.numeric(substr(dates_all, 9, 10))
    out <- data.frame(
      year = Years,
      month = Months,
      day = Days,
      extract = extract_int
    )
  } else {
    out <- data.frame(
      year = Years,
      month = Months,
      extract = extract_int
    )
  }

  if (!downscale) {

    out

  } else {

    ## rescale w/ worldclim climatology

    wc_files <- list.files(worldclim, pattern = ".*\\.bil$")
    wc_param <- strsplit(wc_files[1], "_")[[1]][1]
    spp <- SpatialPoints(
      data.frame(
        x = coords$lon,
        y = coords$lat))

    for (j in month_num) {
      wc <- raster(file.path(worldclim, paste(wc_param, "_", j, ".bil",
                                              sep = "")))
      ## FIXME hier muss der Datentyp richtig geraten werden!
      wc_point <- extract(wc, spp)/scale_fac
      ## when this is NA, we have to look for the next nearest point
      ## that is not NA
      if (is.na(wc_point)) {
        arc30 <-  30/3600
        spp2 <- expand.grid(
          data.frame(
            x = c(coords$lon - arc30, coords$lon + arc30),
            y = c(coords$lat - arc30, coords$lat + arc30)
            ))
        spp2 <- SpatialPoints(spp2)
        wc_points <- extract(wc, spp2)/scale_fac
        wc_point <- wc_points[which(!is.na(wc_points))][1]
      }
      out[out$month == j,]$extract <-
        meth2(out[out$month == j,]$extract, wc_point)
    }

    out

  }
}

##' Approximating grid cell area as trapezoid area
##' @export
gridcellarea <- function(coords, resolution = 0.5) {
  halfres <- resolution/2
  upperleft <- list(
      lon = coords$lon - halfres,
      lat = coords$lat + halfres)
  lowerleft <- list(
      lon = coords$lon - halfres,
      lat = coords$lat - halfres)
  lowerright <- list(
    lon = coords$lon + halfres,
    lat = coords$lat - halfres)
  upperright <- list(
      lon = coords$lon + halfres,
      lat = coords$lat + halfres)
  
  a <- great_circle_dist(lowerleft, lowerright)
  c <- great_circle_dist(upperleft, upperright)
  d <- great_circle_dist(lowerleft, upperleft)
  b <- abs(a - c)/2
  h <- sqrt(d^2 - b^2)
  0.5 * (a + c) * h
}
