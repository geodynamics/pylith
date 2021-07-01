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

#include <portinfo>

#include "TestDataWriterHDF5ExtMesh.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

namespace pylith {
    namespace meshio {

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtMesh_Tri : public TestDataWriterHDF5ExtMesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtMesh_Tri, TestDataWriterHDF5ExtMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtMesh::setUp();
                _data = new TestDataWriterHDF5ExtMesh_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "tri3.h5";
                _data->vertexFilename = "tri3_vertex.h5";
                _data->cellFilename = "tri3_cell.h5";

                TestDataWriterMesh::_setDataTri();
                TestDataWriterMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5ExtMesh_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtMesh_Tri);

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtMesh_Quad : public TestDataWriterHDF5ExtMesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtMesh_Quad, TestDataWriterHDF5ExtMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtMesh::setUp();
                _data = new TestDataWriterHDF5ExtMesh_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "quad4.h5";
                _data->vertexFilename = "quad4_vertex.h5";
                _data->cellFilename = "quad4_cell.h5";

                TestDataWriterMesh::_setDataQuad();
                TestDataWriterMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5ExtMesh_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtMesh_Quad);

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtMesh_Tet : public TestDataWriterHDF5ExtMesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtMesh_Tet, TestDataWriterHDF5ExtMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtMesh::setUp();
                _data = new TestDataWriterHDF5ExtMesh_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "tet4.h5";
                _data->vertexFilename = "tet4_vertex.h5";
                _data->cellFilename = "tet4_cell.h5";

                TestDataWriterMesh::_setDataTet();
                TestDataWriterMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5ExtMesh_Tet
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtMesh_Tet);

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtMesh_Hex : public TestDataWriterHDF5ExtMesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtMesh_Hex, TestDataWriterHDF5ExtMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtMesh::setUp();
                _data = new TestDataWriterHDF5ExtMesh_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "hex8.h5";
                _data->vertexFilename = "hex8_vertex.h5";
                _data->cellFilename = "hex8_cell.h5";

                TestDataWriterMesh::_setDataHex();
                TestDataWriterMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5ExtMesh_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtMesh_Hex);

    } // meshio
} // pylith


// End of file
