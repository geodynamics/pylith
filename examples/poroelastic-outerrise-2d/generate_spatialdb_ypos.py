#!/usr/bin/env nemesis

import numpy

import simpledb_writer

X_WEST = -0.1e+3
X_EAST = 150.1e+3
dx = 0.1e+3 

x_values = numpy.arange(X_WEST, X_EAST + dx, dx)
y_values = numpy.zeros(len(x_values))
x_flexure = x_values

model_run_time = 300e+3
scale_for_flexure = 5e+16

velocity_y = -x_values**4 / scale_for_flexure / model_run_time
velocity_x = numpy.zeros(len(velocity_y))
start_times = numpy.zeros(len(velocity_y))

vertices = numpy.array([x_values, y_values]).T

field_names = ["rate_amplitude_x", "rate_amplitude_y", "rate_start_time"]
field_units = ["m/year", "m/year", "year"]
field_data = numpy.array([[velocity_x], [velocity_y], [start_times]])

filename = 'top_velocity_boundary.spatialdb'

simpledb_writer.write_spatialdb_file(vertices, field_names, field_units, field_data, filename, dimension=2)

rho_water = 1000
g = 9.81
pressure = x_values**4 / scale_for_flexure * g * rho_water
field_names = ["initial_amplitude"]
field_units = ["Pa"]
field_data = numpy.array([[pressure]])

filename = 'surface_fluid_pressure.spatialdb'

simpledb_writer.write_spatialdb_file(vertices, field_names, field_units, field_data, filename, dimension=2)

# End of file
