#!/usr/bin/env python3
"""Python script to extract elevation data in region around Kilamonjaro using GDAL Python API."""

import pathlib

from osgeo import gdal

FILENAME_DEMS_CSV = "30n150w_20101117_gmted_med300.tif,50n150w_20101117_gmted_med300.tif,30n120w_20101117_gmted_med300.tif,50n120w_20101117_gmted_med300.tif"
FILENAME_DEMSUB = "topography.nc"
BBOX = (-130.195718069, 51.5863472565, -115.821471093, 39.0650222389)
RESOLUTION = 0.2 # degrees

FILENAME_DEMS = [filename.strip() for filename in FILENAME_DEMS_CSV.split(",")]

vrt_path = "/vsimem/merged.vrt"
gdal.BuildVRT(vrt_path, FILENAME_DEMS)

gdal.UseExceptions()
options = gdal.TranslateOptions(
        projWin=BBOX,
        xRes=RESOLUTION,
        yRes=RESOLUTION,
        creationOptions=["COMPRESS=DEFLATE"],
        resampleAlg="cubic")
gdal.Translate(FILENAME_DEMSUB, vrt_path, options=options)

gdal.Unlink(vrt_path)
