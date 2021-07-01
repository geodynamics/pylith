#!/usr/bin/env nemesis
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @brief Generate HDF5 data files for eqinfo tests.

import h5py
import numpy

# ----------------------------------------------------------------------


class TestData(object):
    """Abstract base class for test data.
    """def __init__(self):
        self.vertices = None
        self.cells = None
        self.slip = None
        self.time = None
        self.filename = None
        return

    def write(self):
        h5 = h5py.File(self.filename, "w", driver='sec2')
        h5.create_dataset('geometry/vertices', data=self.vertices)
        h5.create_dataset('topology/cells', data=self.cells)
        h5.create_dataset('time', data=self.time)
        slip = h5.create_dataset('vertex_fields/slip', data=self.slip)
        slip.attrs['vector_field_type'] = 'vector'
        h5.close()
        return


# ----------------------------------------------------------------------
class DataLinea(TestData):
    """Test data for fault one with line cells.
    """

    def __init__(self):
        TestData.__init__(self)
        self.vertices = numpy.array([[-1.0, 0.0],
                                     [0.0, 0.0],
                                     [0.0, +1.5]],
                                    dtype=numpy.float64)
        self.cells = numpy.array([[0, 1],
                                  [1, 2]],
                                 dtype=numpy.int32)
        self.slip = numpy.array([
            # t=0
            [[0.4, 0.0],
             [1.0, 0.0],
             [0.8, 0.0]],
            # t=1
            [[0.0, 0.0],
             [0.0, 0.0],
             [0.8, 0.0]]],
            dtype=numpy.float64)
        self.time = numpy.array([0.0, 1.0])
        self.filename = "line_one.h5"


# ----------------------------------------------------------------------
class DataLineb(TestData):
    """Test data for fault two with line cells.
    """

    def __init__(self):
        TestData.__init__(self)
        self.vertices = numpy.array([[-1.5, 0.0],
                                     [0.0, 1.0],
                                     [+1.0, 0.0]],
                                    dtype=numpy.float64)
        self.cells = numpy.array([[1, 2],
                                  [0, 1]],
                                 dtype=numpy.int32)
        self.slip = numpy.array([
            # t=0
            [[0.6, 0.0],
             [1.2, 0.0],
             [0.2, 0.0]],
            # t=1
            [[0.6, 0.0],
             [0.0, 0.0],
             [0.0, 0.0]]],
            dtype=numpy.float64)
        self.time = numpy.array([0.0, 1.0])
        self.filename = "line_two.h5"


# ----------------------------------------------------------------------
class DataTri3a(TestData):
    """Test data for fault one with tri cells.
    """

    def __init__(self):
        TestData.__init__(self)
        self.vertices = numpy.array([[-1.5,  0.0,  0.0],
                                     [0.0,  0.0,  0.0],
                                     [0.0, +2.0,  0.0],
                                     [0.0,  0.0, +2.0]],
                                    dtype=numpy.float64)
        self.cells = numpy.array([[0, 1, 3],
                                  [1, 2, 3]],
                                 dtype=numpy.int32)
        self.slip = numpy.array([
            # t=0
            [[0.0, 0.5, 0.0],
             [0.6, 0.4, 0.0],
             [0.9, 0.2, 0.0],
             [0.0, 0.6, 0.0]],
            # t=1
            [[0.0, 0.3, 0.0],
             [0.0, 0.0, 0.0],
             [0.6, 0.0, 0.0],
             [0.0, 0.0, 0.0]]],
            dtype=numpy.float64)
        self.time = numpy.array([0.0, 1.0])
        self.filename = "tri_one.h5"


# ----------------------------------------------------------------------
class DataTri3b(TestData):
    """Test data for fault two with line cells.
    """

    def __init__(self):
        TestData.__init__(self)
        self.vertices = numpy.array([[-1.5,  0.0,  0.0],
                                     [0.0,  0.0,  0.0],
                                     [0.0,  0.0, +2.0],
                                     [0.0,  2.0,  0.0]],
                                    dtype=numpy.float64)
        self.cells = numpy.array([[0, 1, 3],
                                  [1, 2, 3]],
                                 dtype=numpy.int32)
        self.slip = numpy.array([
            # t=0
            [[0.3, 0.0, 0.0],
             [0.0, 0.0, 0.0],
             [0.0, 0.0, 0.0],
             [0.0, 0.0, 0.0]],
            # t=1
            [[0.0, 0.0, 0.0],
             [0.0, 0.0, 0.0],
             [0.0, 0.0, 0.0],
             [0.0, 0.0, 0.0]]],
            dtype=numpy.float64)
        self.time = numpy.array([0.0, 1.0])
        self.filename = "tri_two.h5"


# ----------------------------------------------------------------------
class DataQuad4a(TestData):
    """Test data for fault one with quad cells.
    """

    def __init__(self):
        TestData.__init__(self)
        self.vertices = numpy.array([[-1.5,  0.0,  0.0],  # 0
                                     [0.0,  0.0,  0.0],  # 1
                                     [0.0, +1.2,  0.0],  # 2
                                     [-1.5,  0.0, +2.0],  # 3
                                     [0.0,  0.0, +1.5],  # 4
                                     [0.0, +1.5, +1.5]],  # 5
                                    dtype=numpy.float64)
        self.cells = numpy.array([[0, 1, 4, 3],
                                  [1, 2, 5, 4]],
                                 dtype=numpy.int32)
        self.slip = numpy.array([
            # t=0
            [[0.1, 1.3, 0.0],
             [0.2, 1.2, 0.0],
             [0.3, 1.1, 0.0],
             [0.4, 1.4, 0.0],
             [0.5, 1.5, 0.0],
             [0.6, 1.6, 0.0]],
            # t=1
            [[0.2, 0.4, 0.0],  # 0
             [0.6, 0.8, 0.0],  # 1
             [1.0, 1.2, 0.0],  # 2
             [1.4, 1.6, 0.0],  # 3
             [1.8, 2.0, 0.0],  # 4
             [2.2, 2.4, 0.0]]],  # 5
            dtype=numpy.float64)
        self.time = numpy.array([0.0, 5.0])
        self.filename = "quad_one.h5"


# ----------------------------------------------------------------------
class DataQuad4b(TestData):
    """Test data for fault two with quad cells.
    """

    def __init__(self):
        TestData.__init__(self)
        self.vertices = numpy.array([[-1.5,  0.0,  0.0],
                                     [0.0,  0.0,  0.0],
                                     [+1.2,  0.0,  0.0],
                                     [-1.5,  0.0, +2.0],
                                     [0.0,  0.0, +1.5],
                                     [+1.5,  0.0, +1.5]],
                                    dtype=numpy.float64)
        self.cells = numpy.array([[0, 1, 4, 3],
                                  [1, 2, 5, 4]],
                                 dtype=numpy.int32)
        self.slip = numpy.array([
            # t=0
            [[0.1, 1.3, 0.0],
             [0.2, 1.2, 0.0],
             [0.3, 1.1, 0.0],
             [0.4, 1.4, 0.0],
             [0.5, 1.5, 0.0],
             [0.6, 1.6, 0.0]],
            # t=1
            [[-0.2, -0.4, 0.0],
             [-0.6, -0.8, 0.0],
             [-1.0, -1.2, 0.0],
             [-1.4, -1.6, 0.0],
             [-1.8, -2.0, 0.0],
             [-2.2, -2.4, 0.0]]],
            dtype=numpy.float64)
        self.time = numpy.array([0.0, 5.0])
        self.filename = "quad_two.h5"


# ======================================================================
# Line mesh
data = DataLinea()
data.write()

data = DataLineb()
data.write()

# Tri mesh
data = DataTri3a()
data.write()

data = DataTri3b()
data.write()

# Quad mesh
data = DataQuad4a()
data.write()

data = DataQuad4b()
data.write()

# End of file
