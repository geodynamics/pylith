import numpy as np
import spatialdb_writer

X_WEST = -0.1e+3
X_EAST = 150.1e+3
dx = 0.1e+3 

x_values = np.arange(X_WEST, X_EAST + dx, dx)
y_values = np.zeros(len(x_values))
x_flexure = x_values

model_run_time = 300e3
scale_for_flexure = 5e16

velocity_y = -x_values**4 / scale_for_flexure / model_run_time
velocity_x = np.zeros(len(velocity_y))
start_times = np.zeros(len(velocity_y))

vertices = np.array([x_values, y_values]).T

field_names = ["rate_amplitude_x", "rate_amplitude_y", "rate_start_time"]

field_units = ["m/year", "m/year", "year"]

field_data = np.array([[velocity_x], [velocity_y], [start_times]])

filename = 'top_velocity_boundary.spatialdb'

spatialdb_writer.write_spatialdb_file(vertices, field_names, field_units, field_data, filename, dimension=2)

rho_water = 1000
g = 9.81
pressure = x_values**4 / scale_for_flexure * g * rho_water
field_names = ["initial_amplitude"]
field_units = ["Pa"]
field_data = np.array([[pressure]])

filename = 'surface_fluid_pressure.spatialdb'

spatialdb_writer.write_spatialdb_file(vertices, field_names, field_units, field_data, filename, dimension=2)

# End file
