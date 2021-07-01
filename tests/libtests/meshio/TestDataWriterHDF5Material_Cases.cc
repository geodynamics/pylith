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

#include "TestDataWriterHDF5Material.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

namespace pylith {
    namespace meshio {

        // --------------------------------------------------------------
        class TestDataWriterHDF5Material_Tri : public TestDataWriterHDF5Material {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5Material_Tri, TestDataWriterHDF5Material);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5Material::setUp();
                _data = new TestDataWriterHDF5Material_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "tri3_mat.h5";
                _data->vertexFilename = "tri3_mat_vertex.h5";
                _data->cellFilename = "tri3_mat_cell.h5";

                TestDataWriterMaterial::_setDataTri();
                TestDataWriterMaterial::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5Material_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5Material_Tri);

        // --------------------------------------------------------------
        class TestDataWriterHDF5Material_Quad : public TestDataWriterHDF5Material {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5Material_Quad, TestDataWriterHDF5Material);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5Material::setUp();
                _data = new TestDataWriterHDF5Material_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "quad4_mat.h5";
                _data->vertexFilename = "quad4_mat_vertex.h5";
                _data->cellFilename = "quad4_mat_cell.h5";

                TestDataWriterMaterial::_setDataQuad();
                TestDataWriterMaterial::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5Material_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5Material_Quad);

        // --------------------------------------------------------------
        class TestDataWriterHDF5Material_Tet : public TestDataWriterHDF5Material {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5Material_Tet, TestDataWriterHDF5Material);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5Material::setUp();
                _data = new TestDataWriterHDF5Material_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "tet4_mat.h5";
                _data->vertexFilename = "tet4_mat_vertex.h5";
                _data->cellFilename = "tet4_mat_cell.h5";

                TestDataWriterMaterial::_setDataTet();
                TestDataWriterMaterial::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5Material_Tet
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5Material_Tet);

        // --------------------------------------------------------------
        class TestDataWriterHDF5Material_Hex : public TestDataWriterHDF5Material {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5Material_Hex, TestDataWriterHDF5Material);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5Material::setUp();
                _data = new TestDataWriterHDF5Material_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "hex8_mat.h5";
                _data->vertexFilename = "hex8_mat_vertex.h5";
                _data->cellFilename = "hex8_mat_cell.h5";

                TestDataWriterMaterial::_setDataHex();
                TestDataWriterMaterial::_initialize();

                PYLITH_METHOD_END;
            } // setUp
        }; // class TestDataWriterHDF5Material_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5Material_Hex);

    } // meshio
} // pylith


// End of file
