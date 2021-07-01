// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file tests/libtests/meshio/TestDataWriter.hh
 *
 * @brief C++ TestDataWriter object
 *
 * C++ unit testing for DataWriter..
 */

#if !defined(pylith_meshio_testdatawriter_hh)
#define pylith_meshio_testdatawriter_hh

#include "pylith/topology/FieldBase.hh" // USES FieldBase::VectorFieldEnum
#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
    namespace meshio {
        class TestDataWriter_Data;
    } // meshio
} // pylith

class pylith::meshio::TestDataWriter_Data {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestDataWriter_Data(void);

    /// Destructor
    ~TestDataWriter_Data(void);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    const char* meshFilename; ///< Name of mesh input file.
    const char* faultLabel; ///< Name of group of vertices for fault.
    int faultId; ///< Material identifier for fault.
    int spaceDim; ///< Spatial dimension.
    PylithReal lengthScale; ///< Length scale for mesh.

    PylithScalar time; ///< Time for fields.
    const char* timeFormat; ///< Format for time stamp.

    PylithInt vertexNumPoints; ///< Number of points in vertex field.
    pylith::topology::FieldBase::Discretization vertexDiscretization; ///< Discretization for vertex fields.
    PylithScalar* vertexValues; ///< Values for vertex field (vertexNumPoints, vertexNumDOF).
    PylithInt vertexNumDOF; ///< Number of values per vertex.

    PylithInt cellNumPoints; ///< Number of points in cell field.
    pylith::topology::FieldBase::Discretization cellDiscretization; ///< Discretization for cell fields.
    PylithScalar* cellValues; ///< Values for cell field (cellNumPoints, cellNumDOF).
    PylithInt cellNumDOF; ///< Number of values per cell.

}; // class TestDataWriter_Data

#endif // pylith_meshio_testdatawriter_hh

// End of file
