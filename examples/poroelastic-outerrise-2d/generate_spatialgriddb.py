import numpy as np
import spatialgrid_writer

X_WEST = -0.1e+3
X_EAST = 150.1e+3
dx = 0.1e+3

Y_BOT = -30.1e+3
Y_TOP = 0.1e+3
dy = 0.1e+3

x_array = np.arange(X_WEST, X_EAST + dx, dx)
y_array = np.arange(Y_BOT, Y_TOP + dy, dy)
z_array = np.array([0.0])
X_matrix, Y_matrix = np.meshgrid(x_array, y_array)
X_flat, Y_flat = X_matrix.flatten(), Y_matrix.flatten()

vs = np.ones(len(X_flat)) * 3000
vp = np.ones(len(X_flat)) * 5300
solid_density = np.ones(len(X_flat)) * 2950 # Turcotte & Schubert Appendix
fluid_density = np.ones(len(X_flat)) * 1000
fluid_viscosity = np.ones(len(X_flat)) * 1e-3
shear_modulus = np.ones(len(X_flat)) * 30e9 # Turcotte & Schubert Appendix
drained_bulk_modulus = np.ones(len(X_flat)) * 10e9
fluid_bulk_modulus = np.ones(len(X_flat)) * 2.0e9
solid_bulk_modulus = np.ones(len(X_flat)) * 80.0e9 # Turcotte & Schubert Appendix (Calculated)
biot_coefficient = 1 - drained_bulk_modulus / solid_bulk_modulus

phi_0 = 0.1
n = 89.53
m = 0.008
porosity = phi_0 / (1 + m * np.abs(Y_flat) / 1e3)**n


###################################### HERE WE SPECIFY THE BACKGROUND PERMEABILITY ######################################

background_permeability = ((solid_density * -9.81 * (Y_flat - 1))**(-4) * 2.5e12)
total_permeability = np.copy(background_permeability)

###################################### HERE WE SPECIFY THE FAULT ZONE PERMEABILITY ######################################

max_fault_depth = -15.0e+3
fault_dips = np.deg2rad(60)

OUTER_RISE_FAULT_SURFACES_X = np.array([100.0e+3, 110.0e+3, 120.0e+3])
OUTER_RISE_FAULT_SURFACES_Y = np.zeros(len(OUTER_RISE_FAULT_SURFACES_X))

OUTER_RISE_FAULT_BURIED_EDGES_Y = np.ones(len(OUTER_RISE_FAULT_SURFACES_X)) * max_fault_depth
OUTER_RISE_FAULT_BURIED_EDGES_X = OUTER_RISE_FAULT_SURFACES_X - OUTER_RISE_FAULT_BURIED_EDGES_Y / np.tan(fault_dips)

surface_points = np.array([OUTER_RISE_FAULT_SURFACES_X, OUTER_RISE_FAULT_SURFACES_Y]).T
buried_points = np.array([OUTER_RISE_FAULT_BURIED_EDGES_X, OUTER_RISE_FAULT_BURIED_EDGES_Y]).T

fault_perm_factor = 1e3
fault_zone_width = 2e3

for i in range(len(surface_points)):
    distance_from_fault = np.abs( (surface_points[i][0] - buried_points[i][0]) * (buried_points[i][1] - Y_flat) - \
                                (buried_points[i][0] - X_flat) * (surface_points[i][1] - buried_points[i][1]) ) / \
                          np.sqrt( (surface_points[i][0] - buried_points[i][0])**2 + (surface_points[i][1] - buried_points[i][1])**2 )
    
    within_fault_zone_index = np.where( (distance_from_fault < fault_zone_width/2) & (Y_flat >= buried_points[i][1]) )
    outside_fault_zone_index = np.where( (distance_from_fault >= fault_zone_width/2) & (Y_flat < buried_points[i][1]) )

    distance_from_fault[outside_fault_zone_index] = fault_zone_width/2

    total_permeability[within_fault_zone_index] = total_permeability[within_fault_zone_index] * fault_perm_factor - \
                                                  (distance_from_fault[within_fault_zone_index] / (fault_zone_width/2)) * \
                                                   total_permeability[within_fault_zone_index] * fault_perm_factor

min_permeability = 1e-26
max_permeability = 1e-12

background_permeability[np.where(background_permeability < min_permeability)] = min_permeability
background_permeability[np.where(background_permeability > max_permeability)] = max_permeability

total_permeability[np.where(total_permeability < min_permeability)] = min_permeability
total_permeability[np.where(total_permeability > max_permeability)] = max_permeability

vertices = np.array([X_flat, Y_flat]).T

field_names = ["solid_density", "fluid_density", "fluid_viscosity", \
               "porosity", "shear_modulus", "drained_bulk_modulus", "biot_coefficient", 
               "fluid_bulk_modulus", "solid_bulk_modulus", "isotropic_permeability"]

field_units = ["kg/m**3", "kg/m**3", "Pa*s", \
               "none", "Pa", "Pa", "none", \
               "Pa", "Pa", "m**2"]

field_data = np.array([[solid_density], [fluid_density], [fluid_viscosity],\
                       [porosity], [shear_modulus], [drained_bulk_modulus], [biot_coefficient], \
                       [fluid_bulk_modulus], [solid_bulk_modulus], [background_permeability]])

field_data_faults= np.array([[solid_density], [fluid_density], [fluid_viscosity],\
                             [porosity], [shear_modulus], [drained_bulk_modulus], [biot_coefficient], \
                             [fluid_bulk_modulus], [solid_bulk_modulus], [total_permeability]])

no_fault_file = 'no_faultzone_permeability.spatialgrid'
fault_file = 'enhanced_faultzone_permeability.spatialgrid'
spatialgrid_writer.write_spatialgrid_file(vertices, x_array, y_array, z_array, field_names, field_units, field_data, no_fault_file, dimension=2)
spatialgrid_writer.write_spatialgrid_file(vertices, x_array, y_array, z_array, field_names, field_units, field_data_faults, fault_file, dimension=2)

# End of file
