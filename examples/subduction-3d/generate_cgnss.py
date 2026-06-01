#!/usr/bin/env nemesis

import numpy
import h5py
import coordsys
from scipy.interpolate import griddata

# Load GNSS station coordinates using NumPy
station_data = numpy.genfromtxt(
    "cgnss_stations_latlon.txt", dtype=None, names=True, encoding="utf-8"
)
station_name = station_data["Site"]
station_latitude = station_data["Latitude"]
station_longitude = station_data["Longitude"]

# Load ground surface from Step 1
h5 = h5py.File("output/step01_axialdisp-groundsurf.h5")
vertices = h5["/geometry/vertices"][:]

station_xy = numpy.zeros((station_name.shape[0], 3), dtype=numpy.float64)
station_xy[:, 0] = station_latitude
station_xy[:, 1] = station_longitude

# Convert from geographic coordinates to mesh coordinates
coordsys.geoToMesh(station_xy)

# Interpolate elevation values from ground surface
station_elevation = griddata(
    vertices[:, 0:2],
    vertices[:, 2],
    (station_xy[:, 0], station_xy[:, 1]),
    method="linear",
)

# Adjust station elevation to slightly below ground surface
station_elevation -= 2.0

# Create array to hold output
data = numpy.empty((station_name.shape[0], 4), dtype=object)
data[:, 0] = station_name
data[:, 1] = station_latitude
data[:, 2] = station_longitude
data[:, 3] = station_elevation

numpy.savetxt(
    "cgnss_stations.txt",
    data,
    fmt="%s %.6f %.6f %.1f",
    header="Site Latitude Longitude Elevation (m)",
)
