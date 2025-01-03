
#### Setup
# Load packages
library(blackmarbler)
library(geodata)
library(sf)
library(raster)
library(ggplot2)

#### Define NASA bearer token
bearer <- "BEARER-TOKEN-HERE"

### ROI
# Define region of interest (roi). The roi must be (1) an sf polygon and (2)
# in the WGS84 (epsg:4326) coordinate reference system. Here, we use the
# getData function to load a polygon of Ghana
roi_sf <- gadm(country = "GHA", level=1, path = tempdir()) |> st_as_sf()
