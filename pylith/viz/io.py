# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
"""Module for reading PyLith output for PyVista."""

import h5py
import numpy
import pyvista

from pylith.viz import core

# Mapping from cell dimension and number of corners to cell type.
CELL_SHAPE = {
    (1, 2): pyvista.CellType.LINE,
    (2, 3): pyvista.CellType.TRIANGLE,
    (2, 4): pyvista.CellType.QUAD,
    (3, 4): pyvista.CellType.TETRA,
    (3, 8): pyvista.CellType.HEXAHEDRON,
}


class PyLithReader:
    """Methods for reading data from PyLith HDF5 output files."""

    @staticmethod
    def read(filename, field_names=None):
        """Read data from PyLith HDF5 output files.
        
        All time steps and fields are read by default.
        """
        h5 = h5py.File(filename)
        grid = PyLithReader._read_grid(h5)
        fields = {}
        fields.update(PyLithReader._read_data_fields(h5, "vertex_fields"))
        fields.update(PyLithReader._read_data_fields(h5, "cell_fields"))
        time = h5["time"][:].squeeze() if h5["time"][:].size > 1 else numpy.array([h5["time"][0]])
        h5.close()
        return core.Mesh(grid=grid, fields=fields, time=time, name=filename)

    @staticmethod
    def _read_grid(h5):
        cells = h5["viz"]["topology"]["cells"]
        cell_dim = cells.attrs["cell_dim"]
        ncorners = cells[:].shape[1]
        cell_type = CELL_SHAPE[(cell_dim, ncorners)]

        vertices = h5["geometry"]["vertices"][:]
        if vertices.shape[1] == 2:
            zeros = numpy.zeros_like(vertices[:,0])
            vertices = numpy.stack((vertices[:,0], vertices[:,1], zeros)).T

        return pyvista.UnstructuredGrid({cell_type: cells[:].astype(int)}, vertices)

    @staticmethod
    def _read_data_fields(h5, group):
        if not group in h5:
            return {}

        data_fields = {}
        basis_order = 1 if group == "vertex_fields" else 0
        for field_name in h5[group]:
            vector_field_type = h5[group][field_name].attrs["vector_field_type"].decode("utf-8")
            if vector_field_type == "vector" and h5[group][field_name].shape[-1] == 2:
                values = h5[group][field_name][:]
                zeros = numpy.zeros_like(values[:,:,0])
                values = numpy.stack((values[:,:,0].T, values[:,:,1].T, zeros.T)).T
            elif vector_field_type == "tensor" and h5[group][field_name].shape[-1] == 4:
                zeros = numpy.zeros_like(values[:,:,0])
                values = numpy.stack((
                    values[:,:,0].T,
                    values[:,:,1].T,
                    values[:,:,2].T,
                    values[:,:,3].T,
                    zeros.T,
                    zeros.T,
                    zeros.T,
                    )).T
            else:
                values = h5[group][field_name][:]

            data_fields[field_name] = core.Field(values=values, field_type=vector_field_type, basis_order=basis_order)
        return data_fields
