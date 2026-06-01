#!/usr/bin/env nemesis
"""Python application to create spatial databases for a synthetic event with time-varying Gaussian slip.

For simplicity we specify the slip on a horizontal plane in latitude and longitude coordinates.
"""

import math
import numpy
from spatialdata.spatialdb.SimpleGridAscii import createWriter
from spatialdata.spatialdb import TimeHistoryIO
import coordsys

FILENAME_SPATIALDB = "fault_slabtop_slowslip.spatialdb"
FILENAME_TIMEDB = "fault_slabtop_slowslip.timedb"

CS_GEO3D = coordsys.cs_geo3D()
LONGITUDE = numpy.arange(-126.4, -123.0, 0.04)
LATITUDE = numpy.arange(44.4, 46.441, 0.04)

RAKE = 110.0
SLIP_CENTER = (45.5, -124.0)
SLIP_MAXRADIUS = 0.7
SLIP_MAX = 8.0
SLIP_STD = 0.3

SLIP_TIME_UNITS = "day"
SLIP_TIME = numpy.array((
    (0.0, 0.0),
    (6.0, 0.2),
    (12.0, 0.5),
    (18.0, 0.8),
    (24.0, 0.9),
    (30.0, 1.0),
    (36.0, 1.0),
))


latitude_2d, longitude_2d = numpy.meshgrid(LATITUDE, LONGITUDE)
elevation_2d = numpy.zeros_like(latitude_2d)
grid = numpy.column_stack((
    latitude_2d.flatten(),
    longitude_2d.flatten(),
    elevation_2d.flatten(),
))

latitude_relative = grid[:, 0] - SLIP_CENTER[0]
longitude_relative = grid[:, 1] - SLIP_CENTER[1]

longitude_normalized = 0.5 * longitude_relative**2 / SLIP_STD**2
latitude_normalized = 0.5 * latitude_relative**2 / SLIP_STD**2

slip = SLIP_MAX * numpy.exp(-(longitude_normalized + latitude_normalized))
distance = numpy.sqrt(longitude_relative**2 + latitude_relative**2)
slip[distance > SLIP_MAXRADIUS] = 0.0

rake_radians = math.radians(RAKE)
left_lateral = slip * math.cos(rake_radians)
reverse_slip = slip * math.sin(rake_radians)
opening_slip = numpy.zeros_like(left_lateral)

slip = numpy.column_stack((left_lateral, reverse_slip, opening_slip))

data = {
    'num-x': LATITUDE.shape[0],
    'num-y': LONGITUDE.shape[0],
    'num-z': 1,
    'points': grid,
    'x': LATITUDE,
    'y': LONGITUDE,
    'z': numpy.array([0.0]),
    'coordsys': CS_GEO3D,
    'data_dim': 2,
    'values': [
        {
            'name': "final_slip_left_lateral",
            'units': "m",
            'data': slip[:, 0]
        },
        {
            'name': "final_slip_reverse",
            'units': "m",
            'data': slip[:, 1]
        },
        {
            'name': "final_slip_opening",
            'units': "m",
            'data': slip[:, 2]
        },
        {
            'name': "initiation_time",
            'units': SLIP_TIME_UNITS,
            'data': numpy.zeros_like(slip[:,0]),
        },
    ]
}
writer = createWriter(FILENAME_SPATIALDB)
writer.write(data)

TimeHistoryIO.write(SLIP_TIME[:,0], SLIP_TIME[:,1], SLIP_TIME_UNITS, FILENAME_TIMEDB)
