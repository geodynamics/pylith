// -*- C++ -*-
//
// -----------------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// -----------------------------------------------------------------------------
//

#include <portinfo>

#include "TestFieldMesh.hh" // Implementation of class methods

// -----------------------------------------------------------------------------
namespace pylith {
    namespace topology {

        // ---------------------------------------------------------------------
        class TestFieldMesh_Quad : public TestFieldMesh {

            CPPUNIT_TEST_SUB_SUITE( TestFieldMesh_Quad, TestFieldMesh );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestFieldMesh::setUp();

                _data->cellDim = 2;
                _data->numVertices = 4;
                _data->numCells = 1;
                _data->numCorners = 4;
                static const int _cells[1*4] = {
                    0, 1, 2, 3,
                };
                _data->cells = const_cast<int*>(_cells);
                static const PylithScalar _coordinates[4*2] = {
                    0.0, 0.0,
                    1.0, 0.0,
                    0.0, 1.0,
                    1.0, 1.0,
                };
                _data->coordinates = const_cast<PylithScalar*>(_coordinates);
            }   // setUp

        };  // TestFieldMesh_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION( TestFieldMesh_Quad );

    }   // topology
}   // pylith


// End of file
