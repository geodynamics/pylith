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

#include "TestDataWriterVTKMaterial.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

namespace pylith {
    namespace meshio {

        // --------------------------------------------------------------
        class TestDataWriterVTKMaterial_Tri : public TestDataWriterVTKMaterial {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKMaterial_Tri, TestDataWriterVTKMaterial);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKMaterial::setUp();
                _data = new TestDataWriterVTKMaterial_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "tri3_mat.vtk";
                _data->vertexFilename = "tri3_mat_vertex.vtk";
                _data->cellFilename = "tri3_mat_cell.vtk";

                TestDataWriterMaterial::_setDataTri();
                TestDataWriterMaterial::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKMaterial_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKMaterial_Tri);

        // --------------------------------------------------------------
        class TestDataWriterVTKMaterial_Quad : public TestDataWriterVTKMaterial {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKMaterial_Quad, TestDataWriterVTKMaterial);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKMaterial::setUp();
                _data = new TestDataWriterVTKMaterial_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "quad4_mat.vtk";
                _data->vertexFilename = "quad4_mat_vertex.vtk";
                _data->cellFilename = "quad4_mat_cell.vtk";

                TestDataWriterMaterial::_setDataQuad();
                TestDataWriterMaterial::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKMaterial_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKMaterial_Quad);

        // --------------------------------------------------------------
        class TestDataWriterVTKMaterial_Tet : public TestDataWriterVTKMaterial {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKMaterial_Tet, TestDataWriterVTKMaterial);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKMaterial::setUp();
                _data = new TestDataWriterVTKMaterial_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "tet4_mat.vtk";
                _data->vertexFilename = "tet4_mat_vertex.vtk";
                _data->cellFilename = "tet4_mat_cell.vtk";

                TestDataWriterMaterial::_setDataTet();
                TestDataWriterMaterial::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKMaterial_Tet
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKMaterial_Tet);

        // --------------------------------------------------------------
        class TestDataWriterVTKMaterial_Hex : public TestDataWriterVTKMaterial {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKMaterial_Hex, TestDataWriterVTKMaterial);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKMaterial::setUp();
                _data = new TestDataWriterVTKMaterial_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "hex8_mat.vtk";
                _data->vertexFilename = "hex8_mat_vertex.vtk";
                _data->cellFilename = "hex8_mat_cell.vtk";

                TestDataWriterMaterial::_setDataHex();
                TestDataWriterMaterial::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterVTKMaterial_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKMaterial_Hex);

    } // meshio
} // pylith


// End of file
