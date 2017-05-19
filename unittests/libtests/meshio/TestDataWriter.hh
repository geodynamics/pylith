// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/meshio/TestDataWriter.hh
 *
 * @brief C++ TestDataWriter object
 *
 * C++ unit testing for DataWriter..
 */

#if !defined(pylith_meshio_testdatawriter_hh)
#define pylith_meshio_testdatawriter_hh

#include "pylith/topology/FieldBase.hh" // USES FieldBase::VectorFieldEnum
#include "pylith/utils/types.hh" // HASA PylithScalar

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriter_Data;
    } // meshio
} // pylith

class pylith::meshio::TestDataWriter_Data {

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Constructor
    TestDataWriter_Data(void);

    /// Destructor
    ~TestDataWriter_Data(void);

    // PUBLIC STRUCTS ///////////////////////////////////////////////////////
public:

    struct FieldInfo {
        const char* name; ///< Name of field
        pylith::topology::FieldBase::VectorFieldEnum vectorFieldType; ///< Type of field.
        const char** components; ///< Names of components.
        int numComponents; ///< Number of field components.
    }; // FieldStruct

    // PUBLIC MEMBERS ///////////////////////////////////////////////////
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
    PylithScalar* vertexScalarValues; ///< Values for scalar vertex field (vertexNumPoints, vertexScalarNumComponents).
    PylithScalar* vertexVectorValues; ///< Values for vector vertex field (vertexNumPoints, vertexVectorNumComponents).
    PylithScalar* vertexTensorValues; ///< Values for tensor vertex field (vertexNumPoints, vertexTensorNumComponents).
    PylithScalar* vertexOtherValues; ///< Values for other vertex field (vertexNumPoints, vertexOtherNumComponents).
    PylithInt vertexScalarNumComponents; ///< Number of components for scalar vertex field.
    PylithInt vertexVectorNumComponents; ///< Number of components for vector vertex field.
    PylithInt vertexTensorNumComponents; ///< Number of components for tensor vertex field.
    PylithInt vertexOtherNumComponents; ///< Number of components for other vertex field.

    PylithInt cellNumPoints; ///< Number of points in cell field.
    pylith::topology::FieldBase::Discretization cellDiscretization; ///< Discretization for cell fields.
    PylithScalar* cellScalarValues; ///< Values for scalar cell field (cellNumPoints, cellScalarNumComponents).
    PylithScalar* cellVectorValues; ///< Values for vector cell field (cellNumPoints, cellVectorNumComponents).
    PylithScalar* cellTensorValues; ///< Values for tensor cell field (cellNumPoints, cellTensorNumComponents).
    PylithScalar* cellOtherValues; ///< Values for other cell field (cellNumPoints, cellOtherNumComponents).
    PylithInt cellScalarNumComponents; ///< Number of components for scalar cell field.
    PylithInt cellVectorNumComponents; ///< Number of components for vector cell field.
    PylithInt cellTensorNumComponents; ///< Number of components for tensor cell field.
    PylithInt cellOtherNumComponents; ///< Number of components for other cell field.


}; // class TestDataWriter_Data


#endif // pylith_meshio_testdatawriter_hh


// End of file
