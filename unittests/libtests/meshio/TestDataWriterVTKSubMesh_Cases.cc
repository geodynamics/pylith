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

#include <portinfo>

#include "TestDataWriterVTKSubMesh.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

namespace pylith {
    namespace meshio {

        // --------------------------------------------------------------
        class TestDataWriterVTKSubMesh_Tri : public TestDataWriterVTKSubMesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKSubMesh_Tri, TestDataWriterVTKSubMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKSubMesh::setUp();
                _data = new TestDataWriterVTKSubMesh_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "tri3_surf.vtk";
                _data->vertexFilename = "tri3_surf_vertex.vtk";
                _data->cellFilename = "tri3_surf_cell.vtk";

                TestDataWriterSubMesh::_setDataTri();
                TestDataWriterSubMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKSubMesh_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKSubMesh_Tri);

        // --------------------------------------------------------------
        class TestDataWriterVTKSubMesh_Quad : public TestDataWriterVTKSubMesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKSubMesh_Quad, TestDataWriterVTKSubMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKSubMesh::setUp();
                _data = new TestDataWriterVTKSubMesh_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "quad4_surf.vtk";
                _data->vertexFilename = "quad4_surf_vertex.vtk";
                _data->cellFilename = "quad4_surf_cell.vtk";

                TestDataWriterSubMesh::_setDataQuad();
                TestDataWriterSubMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKSubMesh_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKSubMesh_Quad);

        // --------------------------------------------------------------
        class TestDataWriterVTKSubMesh_Tet : public TestDataWriterVTKSubMesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKSubMesh_Tet, TestDataWriterVTKSubMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKSubMesh::setUp();
                _data = new TestDataWriterVTKSubMesh_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "tet4_surf.vtk";
                _data->vertexFilename = "tet4_surf_vertex.vtk";
                _data->cellFilename = "tet4_surf_cell.vtk";

                TestDataWriterSubMesh::_setDataTet();
                TestDataWriterSubMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKSubMesh_Tet
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKSubMesh_Tet);

        // --------------------------------------------------------------
        class TestDataWriterVTKSubMesh_Hex : public TestDataWriterVTKSubMesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKSubMesh_Hex, TestDataWriterVTKSubMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKSubMesh::setUp();
                _data = new TestDataWriterVTKSubMesh_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "hex8_surf.vtk";
                _data->vertexFilename = "hex8_surf_vertex.vtk";
                _data->cellFilename = "hex8_surf_cell.vtk";

                TestDataWriterSubMesh::_setDataHex();
                TestDataWriterSubMesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKSubMesh_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKSubMesh_Hex);

    } // meshio
} // pylith


// End of file
