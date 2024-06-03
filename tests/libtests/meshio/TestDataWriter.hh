// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

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

// End of file
