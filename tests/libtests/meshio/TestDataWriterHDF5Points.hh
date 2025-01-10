// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "TestDataWriterHDF5.hh" // ISA TestDataWriterHDF5
#include "TestDataWriterPoints.hh" // ISA TestDataWriterPoints

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5Points;
        class TestDataWriterHDF5Points_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5Points : public TestDataWriterHDF5, public TestDataWriterPoints {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDataWriterHDF5Points(TestDataWriterHDF5Points_Data* data);

    /// Destructor.
    ~TestDataWriterHDF5Points(void);

    /// Test openTimeStep() and closeTimeStep()
    void testOpenClose(void);

    /// Test writeVertexField.
    void testWriteVertexField(void);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get test data.
     *
     * @returns Test data.
     */
    TestDataWriterPoints_Data* _getData(void);

    // PROTECTED MEMBDERS /////////////////////////////////////////////////
protected:

    TestDataWriterHDF5Points_Data* _data; ///< Data for testing.

}; // class TestDataWriterHDF5Points

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5Points_Data : public TestDataWriterHDF5_Data, public TestDataWriterPoints_Data {};

// End of file
