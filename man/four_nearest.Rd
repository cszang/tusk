\name{four_nearest}
\alias{four_nearest}
\title{Get four nearest gridpoints}
\usage{
  four_nearest(coords, netcdf, data_set = "ts320")
}
\arguments{
  \item{coords}{the coordinates to find the nearest
  gridpoints for as list with $lon and $lat in decimal
  degrees}

  \item{netcdf}{the ncdf object holding the relevant data}

  \item{data_set}{the kind of data used (see
  ?supported_sets) for details}
}
\value{
  a list holding four lists, with each holding a lat/lot
}
\description{
  Get four nearest gridpoints for a lon/lat pair.
}
\details{
  The returned list of gridpoint coordinates can be used
  for interpolation
}

