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

#include "TestDataWriterHDF5ExtMaterial.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

namespace pylith {
    namespace meshio {

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtMaterial_Tri : public TestDataWriterHDF5ExtMaterial {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtMaterial_Tri, TestDataWriterHDF5ExtMaterial);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtMaterial::setUp();
                _data = new TestDataWriterHDF5ExtMaterial_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "tri3_mat.h5";
                _data->vertexFilename = "tri3_mat_vertex.h5";
                _data->cellFilename = "tri3_mat_cell.h5";

                TestDataWriterMaterial::_setDataTri();
                TestDataWriterMaterial::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5ExtMaterial_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtMaterial_Tri);

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtMaterial_Quad : public TestDataWriterHDF5ExtMaterial {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtMaterial_Quad, TestDataWriterHDF5ExtMaterial);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtMaterial::setUp();
                _data = new TestDataWriterHDF5ExtMaterial_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "quad4_mat.h5";
                _data->vertexFilename = "quad4_mat_vertex.h5";
                _data->cellFilename = "quad4_mat_cell.h5";

                TestDataWriterMaterial::_setDataQuad();
                TestDataWriterMaterial::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5ExtMaterial_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtMaterial_Quad);

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtMaterial_Tet : public TestDataWriterHDF5ExtMaterial {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtMaterial_Tet, TestDataWriterHDF5ExtMaterial);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtMaterial::setUp();
                _data = new TestDataWriterHDF5ExtMaterial_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "tet4_mat.h5";
                _data->vertexFilename = "tet4_mat_vertex.h5";
                _data->cellFilename = "tet4_mat_cell.h5";

                TestDataWriterMaterial::_setDataTet();
                TestDataWriterMaterial::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5ExtMaterial_Tet
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtMaterial_Tet);

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtMaterial_Hex : public TestDataWriterHDF5ExtMaterial {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtMaterial_Hex, TestDataWriterHDF5ExtMaterial);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtMaterial::setUp();
                _data = new TestDataWriterHDF5ExtMaterial_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "hex8_mat.h5";
                _data->vertexFilename = "hex8_mat_vertex.h5";
                _data->cellFilename = "hex8_mat_cell.h5";

                TestDataWriterMaterial::_setDataHex();
                TestDataWriterMaterial::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5ExtMaterial_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtMaterial_Hex);

    } // meshio
} // pylith


// End of file
