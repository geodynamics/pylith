#!/usr/bin/env nemesis
"""Python script to extract elevation data in region around the Cascadia subduction zone
using the GDAL Python API.

Notes:

1. To use this script, you need to have downloaded the digital elevation files from
   https://earthexplorer.usgs.gov and unpacked the zip files in the `input` directory.

2. This Python script requires the GDAL (https://gdal.org) Python package and
   its underlying C++ library.
"""
import pathlib

from osgeo import gdal

gdal.UseExceptions()

DATA_DIR = pathlib.Path("input")

FILENAMES_IN = (
    DATA_DIR / "30n120w_20101117_gmted_med300.tif",
    DATA_DIR / "30n150w_20101117_gmted_med300.tif",
    DATA_DIR / "50n120w_20101117_gmted_med300.tif",
    DATA_DIR / "50n150w_20101117_gmted_med300.tif",
)
FILENAME_OUT = "topography.nc"
VRT_PATH = "/vsimem/merged.vrt"

BBOX = (-130.2, 51.6, -115.8, 39.0)
RESOLUTION = 0.2  # degrees

gdal.BuildVRT(VRT_PATH, FILENAMES_IN)
options = gdal.TranslateOptions(
    projWin=BBOX,
    xRes=RESOLUTION,
    yRes=RESOLUTION,
    creationOptions=["COMPRESS=DEFLATE"],
    resampleAlg="cubic",
)
gdal.Translate(FILENAME_OUT, VRT_PATH, options=options)
gdal.Unlink(VRT_PATH)
