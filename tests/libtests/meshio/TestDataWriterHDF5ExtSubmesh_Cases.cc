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

#include "TestDataWriterHDF5ExtSubmesh.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

namespace pylith {
    namespace meshio {
        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtSubmesh_Tri : public TestDataWriterHDF5ExtSubmesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtSubmesh_Tri, TestDataWriterHDF5ExtSubmesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtSubmesh::setUp();
                _data = new TestDataWriterHDF5ExtSubmesh_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "tri3_surf.h5";
                _data->vertexFilename = "tri3_surf_vertex.h5";
                _data->cellFilename = "tri3_surf_cell.h5";

                TestDataWriterSubmesh::_setDataTri();
                TestDataWriterSubmesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp

        }; // class TestDataWriterHDF5ExtSubmesh_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtSubmesh_Tri);

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtSubmesh_Quad : public TestDataWriterHDF5ExtSubmesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtSubmesh_Quad, TestDataWriterHDF5ExtSubmesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtSubmesh::setUp();
                _data = new TestDataWriterHDF5ExtSubmesh_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "quad4_surf.h5";
                _data->vertexFilename = "quad4_surf_vertex.h5";
                _data->cellFilename = "quad4_surf_cell.h5";

                TestDataWriterSubmesh::_setDataQuad();
                TestDataWriterSubmesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp

        }; // class TestDataWriterHDF5ExtSubmesh_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtSubmesh_Quad);

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtSubmesh_Tet : public TestDataWriterHDF5ExtSubmesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtSubmesh_Tet, TestDataWriterHDF5ExtSubmesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtSubmesh::setUp();
                _data = new TestDataWriterHDF5ExtSubmesh_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "tet4_surf.h5";
                _data->vertexFilename = "tet4_surf_vertex.h5";
                _data->cellFilename = "tet4_surf_cell.h5";

                TestDataWriterSubmesh::_setDataTet();
                TestDataWriterSubmesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp

        }; // class TestDataWriterHDF5ExtSubmesh_Tet
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtSubmesh_Tet);

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtSubmesh_Hex : public TestDataWriterHDF5ExtSubmesh {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtSubmesh_Hex, TestDataWriterHDF5ExtSubmesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtSubmesh::setUp();
                _data = new TestDataWriterHDF5ExtSubmesh_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "hex8_surf.h5";
                _data->vertexFilename = "hex8_surf_vertex.h5";
                _data->cellFilename = "hex8_surf_cell.h5";

                TestDataWriterSubmesh::_setDataHex();
                TestDataWriterSubmesh::_initialize();

                PYLITH_METHOD_END;
            } // setUp

        }; // class TestDataWriterHDF5ExtSubmesh_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtSubmesh_Hex);

    } // meshio
} // pylith

// End of file
