from spatialdata.spatialdb.SimpleIOAscii import createWriter
from spatialdata.geocoords.CSCart import CSCart
import numpy

def fault_slip(x_values, y_values, maximum_depth_of_slip, minimum_depth_of_slip, specified_slip_points, polynomial_degree):
    '''
    Creates an array for fault slip, which is fed to write_spatialdb_file to create spatialDB files defining slip on a fault
    in the PyLith model. x_values, y_values are the spatial coordinates defining the profile where slip will be specified. 
    maximum_depth_of_slip and minimum_depth_of_slip specify the depth extents over which slip can be nonzero. specified_slip_points
    is an array of points where slip is specified, between specified points the slip is fit using a polynomial with degree set by
    polynomial_degree.
    '''
    slip_array = numpy.zeros(len(x_values))
    polynomial_coefficients = numpy.polyfit(specified_slip_points[:, 0], specified_slip_points[:, 1], deg=polynomial_degree)

    for i in range(len(y_values)):
        if y_values[i] <= minimum_depth_of_slip and y_values[i] >= maximum_depth_of_slip:
            slip_value = 0
            for p in range(len(polynomial_coefficients)):
                slip_value += polynomial_coefficients[p] * y_values[i]**(len(polynomial_coefficients) - 1 - p)
            slip_array[i] = slip_value
        else:
            slip_array[i] = 0.0
    return slip_array

def write_spatialdb_file(vertices, field_names, field_units, field_data, filename, dimension=2):
    '''
    Writes the spatialDB file. vertices is an nxdimension array containing the spatial coordinates of the spatialDB file.
    rupture_times, left_lateral_slips, opening_angles correspond to when the fault will rupture, how much it will slip,
    and how much it will open, respectively. filename is the name of the spatialDB file.
    '''
    cs = CSCart()
    cs._configure()
    cs.setSpaceDim(dimension)

    spatialdb_string = []
    for i in range(len(field_names)):
        spatialdb_string.append({'name':field_names[i], 'units':field_units[i], 'data':field_data[i]})
   #  print(spatialdb_string)
    writer = createWriter(filename)
    writer.write({'points': vertices,
                  'coordsys': cs,
                  'data_dim': 1,
                  'values': spatialdb_string})
