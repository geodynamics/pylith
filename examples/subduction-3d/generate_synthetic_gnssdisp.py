#!/usr/bin/env nemesis
"""Python application to create synthetic data from PyLith points output.
This script is meant to be used with example step07. It reads the
slow slip displacements at synthetic cGNSS stations from example
step06, adds noise to the displacements, and then generates a text file
with the synthetic data.
"""

import pathlib

import numpy
import h5py

time_step = -1  # final time step
point_input_file = pathlib.Path("output") / "step06_slowslip-cgnss_stations.h5"
output_file = "cgnss_synthetic_displacement.txt"

output_header = ",".join(
    ["Station", "X", "Y", "Z", "UEast", "UNorth", "UUp", "SigEast", "SigNorth", "SigUp"]
)
output_format = ",".join(["%s"] + 9 * ["%g"]) + "\n"

# We assume data uncertainties of 0.5 mm horizontal and 1 mm vertical.
sigma_east = 0.0005
sigma_north = 0.0005
sigma_up = 0.001

h5 = h5py.File(point_input_file, "r")
coords = h5["geometry/vertices"][:]
stations = h5["stations"][:].astype(str)
disp_raw = h5["vertex_fields/displacement"][time_step, :, :]
h5.close()

n_stations = coords.shape[0]

disp_noise = disp_raw.copy()
disp_noise[:, 0] += sigma_east * numpy.random.randn(n_stations)
disp_noise[:, 1] += sigma_north * numpy.random.randn(n_stations)
disp_noise[:, 2] += sigma_up * numpy.random.randn(n_stations)

with open(output_file, "w", encoding="utf-8") as fout:
    fout.write(output_header + "\n")

    for i_station in range(n_stations):
        outLine = output_format % (
            stations[i_station],
            coords[i_station, 0],
            coords[i_station, 1],
            coords[i_station, 2],
            disp_noise[i_station, 0],
            disp_noise[i_station, 1],
            disp_noise[i_station, 2],
            sigma_east,
            sigma_north,
            sigma_up,
        )
        fout.write(outLine)
