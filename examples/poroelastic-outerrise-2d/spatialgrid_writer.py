from spatialdata.spatialdb.SimpleGridAscii import createWriter
from spatialdata.geocoords.CSCart import CSCart
import numpy as np



def write_spatialgrid_file(vertices, x_points, y_points, z_points, field_names, field_units, field_data, filename, dimension=2):
    '''
    Writes a SpatialGridDB file. vertices is an n x dimension array containing the spatial coordinates of the SpatialGridDB file.
    x_points, y_points, z_points are the points that discretize the x-direction, y-direction, and z-direction, respectively.
    field_names is an array that contains the names of the fields that the SpatialGridDB file initializes in the PyLith model.
    field_units is an array that contains the units corresponding to the entries of field_names that the SpatialGridDB file 
    initializes in the PyLith model.
    field_data is an array that contains the data corresponding to the entires of field_names that the SpatialGridDB file
    initializes in the PyLith model. 
    filename is the name of the SpatialGridDB file.
    dimension is the dimension of the SpatialGridDB file.
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
                  'x': x_points,
                  'y': y_points,
                  'z': z_points,
                  'coordsys': cs,
                  'data_dim': 2,
                  'values': spatialdb_string})
# End File