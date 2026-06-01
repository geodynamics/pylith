#!/usr/bin/env python3
"""Python script to pre-process bathymetry and topography data in region around the Cascadia subduction zone
using rasterio.

Notes:

1. To use this script, you need to have extracted the digital elevation file from
   https://www.ncei.noaa.gov/products/etopo-global-relief-model and saved the GeoTiff
   file to the `input` directory.

   Bounding box: -128.5, -117.0; 41.5, 49.5

2. This Python script requires the rasterio Python package, which can be installed
   using pip (`pip install rasterio`).
"""

import pathlib

import rasterio
import rasterio.warp
import netCDF4
import numpy

FILENAME_IN = pathlib.Path("input") / "etopo2020_bedrock_local.tiff"
FILENAME_OUT = "etopo2020_bedrock_local.nc"
DECIMATE = 25  # Decimate the resolution to get a coarse grid for this example

with rasterio.open(FILENAME_IN) as src:
    data = src.read(1)  # Read the first band in the GeoTiff file

    # Create affine transformation with decimated size
    crs = src.crs
    height_new = data.shape[0] // DECIMATE + 1
    width_new = data.shape[1] // DECIMATE + 1
    transform_new = rasterio.Affine(
        src.transform.a * DECIMATE,
        src.transform.b,
        src.transform.c,
        src.transform.d,
        src.transform.e * DECIMATE,
        src.transform.f,
    )

    # Sample data to decimated size
    data_new = numpy.empty((height_new, width_new), dtype=data.dtype)
    rasterio.warp.reproject(
        source=data,
        destination=data_new,
        src_transform=src.transform,
        src_crs=crs,
        dst_transform=transform_new,
        dst_crs=crs,
        resampling=rasterio.enums.Resampling.bilinear,  # or Resampling.nearest for integer values
    )

    # Get coordinates of points using transformation
    latitude = transform_new[5] + numpy.arange(height_new) * transform_new[4]
    longitude = transform_new[2] + numpy.arange(width_new) * transform_new[0]

# Write NetCDF file
with netCDF4.Dataset(FILENAME_OUT, "w", format="NETCDF4") as nc:
    # Create dimensions
    nc.createDimension("lat", height_new)
    nc.createDimension("lon", width_new)

    # Create variables
    crs_var = nc.createVariable("crs", "S1")
    latitude_var = nc.createVariable("lat", numpy.float32, ("lat",))
    longitude_var = nc.createVariable("lon", numpy.float32, ("lon",))
    data_var = nc.createVariable("data", numpy.float32, ("lat", "lon"))

    # Populate variable data, including attributes
    crs_var[:] = ""
    crs_var.grid_mapping_name = "latitude_longitude"
    crs_var.long_name = "CRS definition"
    crs_var.longitude_of_prime_meridian = 0.0
    crs_var.semi_major_axis = 6378137.0
    crs_var.inverse_flattening = 298.257223563
    crs_var.spatial_ref = crs.to_wkt()
    crs_var.crs_wkt = crs.to_wkt()
    crs_var.geoTransform = f"{transform_new[2]} {transform_new[0]} {transform_new[1]} {transform_new[5]} {transform_new[3]} {transform_new[4]}"

    latitude_var[:] = latitude
    latitude_var.standard_name = "latitude"
    latitude_var.long_name = "latitude"
    latitude_var.units = "degrees_north"

    longitude_var[:] = longitude
    longitude_var.standard_name = "longitude"
    longitude_var.long_name = "longitude"
    longitude_var.units = "degrees_east"

    data_var[:, :] = data_new
    data_var.units = "meter"
    data_var.long_name = "Elevation of bedrock."
