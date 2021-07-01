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

#include "TestDataWriterHDF5Mesh.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

namespace pylith {
    namespace meshio {

        // --------------------------------------------------------------
        class TestDataWriterHDF5Mesh_Tri : public TestDataWriterHDF5Mesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5Mesh_Tri, TestDataWriterHDF5Mesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5Mesh::setUp();
                _data = new TestDataWriterHDF5Mesh_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "tri3.h5";
                _data->vertexFilename = "tri3_vertex.h5";
                _data->cellFilename = "tri3_cell.h5";

                TestDataWriterMesh::_setDataTri();
                TestDataWriterMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5Mesh_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5Mesh_Tri);

        // --------------------------------------------------------------
        class TestDataWriterHDF5Mesh_Quad : public TestDataWriterHDF5Mesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5Mesh_Quad, TestDataWriterHDF5Mesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5Mesh::setUp();
                _data = new TestDataWriterHDF5Mesh_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "quad4.h5";
                _data->vertexFilename = "quad4_vertex.h5";
                _data->cellFilename = "quad4_cell.h5";

                TestDataWriterMesh::_setDataQuad();
                TestDataWriterMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5Mesh_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5Mesh_Quad);

        // --------------------------------------------------------------
        class TestDataWriterHDF5Mesh_Tet : public TestDataWriterHDF5Mesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5Mesh_Tet, TestDataWriterHDF5Mesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5Mesh::setUp();
                _data = new TestDataWriterHDF5Mesh_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "tet4.h5";
                _data->vertexFilename = "tet4_vertex.h5";
                _data->cellFilename = "tet4_cell.h5";

                TestDataWriterMesh::_setDataTet();
                TestDataWriterMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5Mesh_Tet
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5Mesh_Tet);

        // --------------------------------------------------------------
        class TestDataWriterHDF5Mesh_Hex : public TestDataWriterHDF5Mesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5Mesh_Hex, TestDataWriterHDF5Mesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5Mesh::setUp();
                _data = new TestDataWriterHDF5Mesh_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "hex8.h5";
                _data->vertexFilename = "hex8_vertex.h5";
                _data->cellFilename = "hex8_cell.h5";

                TestDataWriterMesh::_setDataHex();
                TestDataWriterMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5Mesh_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5Mesh_Hex);

    } // meshio
} // pylith


// End of file
