#!/usr/bin/env python3
"""Python script to extract elevation data in region around Kilamonjaro using GDAL Python API."""

import pathlib

from osgeo import gdal

FILENAME_DEM = "30n150w_20101117_gmted_med300.tif"
FILENAME_DEMSUB = "topography.nc"
BBOX = (-128.195718069, 51.5863472565, -120.821471093, 39.0650222389)
RESOLUTION = 0.1 # degrees


gdal.UseExceptions()
options = gdal.TranslateOptions(
        projWin=BBOX,
        xRes=RESOLUTION,
        yRes=RESOLUTION,
        creationOptions=["COMPRESS=DEFLATE"],
        resampleAlg="cubic")
gdal.Translate(FILENAME_DEMSUB, str(FILENAME_DEM), options=options)
