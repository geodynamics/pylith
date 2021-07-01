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

#include "TestDataWriterVTKMesh.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

namespace pylith {
    namespace meshio {

        // --------------------------------------------------------------
        class TestDataWriterVTKMesh_Tri : public TestDataWriterVTKMesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKMesh_Tri, TestDataWriterVTKMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKMesh::setUp();
                _data = new TestDataWriterVTKMesh_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "tri3.vtk";
                _data->vertexFilename = "tri3_vertex.vtk";
                _data->cellFilename = "tri3_cell.vtk";

                TestDataWriterMesh::_setDataTri();
                TestDataWriterMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKMesh_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKMesh_Tri);

        // --------------------------------------------------------------
        class TestDataWriterVTKMesh_Quad : public TestDataWriterVTKMesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKMesh_Quad, TestDataWriterVTKMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKMesh::setUp();
                _data = new TestDataWriterVTKMesh_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "quad4.vtk";
                _data->vertexFilename = "quad4_vertex.vtk";
                _data->cellFilename = "quad4_cell.vtk";

                TestDataWriterMesh::_setDataQuad();
                TestDataWriterMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKMesh_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKMesh_Quad);

        // --------------------------------------------------------------
        class TestDataWriterVTKMesh_Tet : public TestDataWriterVTKMesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKMesh_Tet, TestDataWriterVTKMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKMesh::setUp();
                _data = new TestDataWriterVTKMesh_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "tet4.vtk";
                _data->vertexFilename = "tet4_vertex.vtk";
                _data->cellFilename = "tet4_cell.vtk";

                TestDataWriterMesh::_setDataTet();
                TestDataWriterMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKMesh_Tet
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKMesh_Tet);

        // --------------------------------------------------------------
        class TestDataWriterVTKMesh_Hex : public TestDataWriterVTKMesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKMesh_Hex, TestDataWriterVTKMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKMesh::setUp();
                _data = new TestDataWriterVTKMesh_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "hex8.vtk";
                _data->vertexFilename = "hex8_vertex.vtk";
                _data->cellFilename = "hex8_cell.vtk";

                TestDataWriterMesh::_setDataHex();
                TestDataWriterMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKMesh_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKMesh_Hex);

    } // meshio
} // pylith


// End of file
