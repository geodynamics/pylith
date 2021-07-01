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

#include "TestDataWriterVTKSubmesh.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

namespace pylith {
    namespace meshio {

        // --------------------------------------------------------------
        class TestDataWriterVTKSubmesh_Tri : public TestDataWriterVTKSubmesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKSubmesh_Tri, TestDataWriterVTKSubmesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKSubmesh::setUp();
                _data = new TestDataWriterVTKSubmesh_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "tri3_surf.vtk";
                _data->vertexFilename = "tri3_surf_vertex.vtk";
                _data->cellFilename = "tri3_surf_cell.vtk";

                TestDataWriterSubmesh::_setDataTri();
                TestDataWriterSubmesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKSubmesh_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKSubmesh_Tri);

        // --------------------------------------------------------------
        class TestDataWriterVTKSubmesh_Quad : public TestDataWriterVTKSubmesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKSubmesh_Quad, TestDataWriterVTKSubmesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKSubmesh::setUp();
                _data = new TestDataWriterVTKSubmesh_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "quad4_surf.vtk";
                _data->vertexFilename = "quad4_surf_vertex.vtk";
                _data->cellFilename = "quad4_surf_cell.vtk";

                TestDataWriterSubmesh::_setDataQuad();
                TestDataWriterSubmesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKSubmesh_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKSubmesh_Quad);

        // --------------------------------------------------------------
        class TestDataWriterVTKSubmesh_Tet : public TestDataWriterVTKSubmesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKSubmesh_Tet, TestDataWriterVTKSubmesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKSubmesh::setUp();
                _data = new TestDataWriterVTKSubmesh_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "tet4_surf.vtk";
                _data->vertexFilename = "tet4_surf_vertex.vtk";
                _data->cellFilename = "tet4_surf_cell.vtk";

                TestDataWriterSubmesh::_setDataTet();
                TestDataWriterSubmesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKSubmesh_Tet
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKSubmesh_Tet);

        // --------------------------------------------------------------
        class TestDataWriterVTKSubmesh_Hex : public TestDataWriterVTKSubmesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKSubmesh_Hex, TestDataWriterVTKSubmesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKSubmesh::setUp();
                _data = new TestDataWriterVTKSubmesh_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "hex8_surf.vtk";
                _data->vertexFilename = "hex8_surf_vertex.vtk";
                _data->cellFilename = "hex8_surf_cell.vtk";

                TestDataWriterSubmesh::_setDataHex();
                TestDataWriterSubmesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKSubmesh_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKSubmesh_Hex);

    } // meshio
} // pylith


// End of file
