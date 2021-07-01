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

#include "TestDataWriterHDF5ExtPoints.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

namespace pylith {
    namespace meshio {
        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtPoints_Tri : public TestDataWriterHDF5ExtPoints {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtPoints_Tri, TestDataWriterHDF5ExtPoints);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtPoints::setUp();
                _data = new TestDataWriterHDF5ExtPoints_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "tri3_points.h5";
                _data->vertexFilename = "tri3_points_vertex.h5";

                TestDataWriterPoints::_setDataTri();
                TestDataWriterPoints::_initialize();

                PYLITH_METHOD_END;
            } // setUp

        }; // class TestDataWriterHDF5ExtPoints_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtPoints_Tri);

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtPoints_Quad : public TestDataWriterHDF5ExtPoints {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtPoints_Quad, TestDataWriterHDF5ExtPoints);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtPoints::setUp();
                _data = new TestDataWriterHDF5ExtPoints_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "quad4_points.h5";
                _data->vertexFilename = "quad4_points_vertex.h5";

                TestDataWriterPoints::_setDataQuad();
                TestDataWriterPoints::_initialize();

                PYLITH_METHOD_END;
            } // setUp

        }; // class TestDataWriterHDF5ExtPoints_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtPoints_Quad);

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtPoints_Tet : public TestDataWriterHDF5ExtPoints {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtPoints_Tet, TestDataWriterHDF5ExtPoints);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtPoints::setUp();
                _data = new TestDataWriterHDF5ExtPoints_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "tet4_points.h5";
                _data->vertexFilename = "tet4_points_vertex.h5";

                TestDataWriterPoints::_setDataTet();
                TestDataWriterPoints::_initialize();

                PYLITH_METHOD_END;
            } // setUp

        }; // class TestDataWriterHDF5ExtPoints_Tet
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtPoints_Tet);

        // --------------------------------------------------------------
        class TestDataWriterHDF5ExtPoints_Hex : public TestDataWriterHDF5ExtPoints {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterHDF5ExtPoints_Hex, TestDataWriterHDF5ExtPoints);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterHDF5ExtPoints::setUp();
                _data = new TestDataWriterHDF5ExtPoints_Data();CPPUNIT_ASSERT(_data);

                _data->opencloseFilename = "hex8_points.h5";
                _data->vertexFilename = "hex8_points_vertex.h5";

                TestDataWriterPoints::_setDataHex();
                TestDataWriterPoints::_initialize();

                PYLITH_METHOD_END;
            } // setUp

        }; // class TestDataWriterHDF5ExtPoints_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterHDF5ExtPoints_Hex);

    } // meshio
} // pylith

// End of file
