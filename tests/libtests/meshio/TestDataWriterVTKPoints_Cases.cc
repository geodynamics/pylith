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

#include "TestDataWriterVTKPoints.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

namespace pylith {
    namespace meshio {
        // --------------------------------------------------------------
        class TestDataWriterVTKPoints_Tri : public TestDataWriterVTKPoints {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKPoints_Tri, TestDataWriterVTKPoints);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKPoints::setUp();
                _data = new TestDataWriterVTKPoints_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "tri3_points.vtk";
                _data->vertexFilename = "tri3_points_vertex.vtk";

                TestDataWriterPoints::_setDataTri();
                TestDataWriterPoints::_initialize();

                PYLITH_METHOD_END;
            } // setUp

        }; // class TestDataWriterVTKPoints_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKPoints_Tri);

        // --------------------------------------------------------------
        class TestDataWriterVTKPoints_Quad : public TestDataWriterVTKPoints {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKPoints_Quad, TestDataWriterVTKPoints);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKPoints::setUp();
                _data = new TestDataWriterVTKPoints_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "quad4_points.vtk";
                _data->vertexFilename = "quad4_points_vertex.vtk";

                TestDataWriterPoints::_setDataQuad();
                TestDataWriterPoints::_initialize();

                PYLITH_METHOD_END;
            } // setUp

        }; // class TestDataWriterVTKPoints_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKPoints_Quad);

        // --------------------------------------------------------------
        class TestDataWriterVTKPoints_Tet : public TestDataWriterVTKPoints {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKPoints_Tet, TestDataWriterVTKPoints);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKPoints::setUp();
                _data = new TestDataWriterVTKPoints_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "tet4_points.vtk";
                _data->vertexFilename = "tet4_points_vertex.vtk";

                TestDataWriterPoints::_setDataTet();
                TestDataWriterPoints::_initialize();

                PYLITH_METHOD_END;
            } // setUp

        }; // class TestDataWriterVTKPoints_Tet
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKPoints_Tet);

        // --------------------------------------------------------------
        class TestDataWriterVTKPoints_Hex : public TestDataWriterVTKPoints {
            CPPUNIT_TEST_SUB_SUITE(TestDataWriterVTKPoints_Hex, TestDataWriterVTKPoints);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                PYLITH_METHOD_BEGIN;

                TestDataWriterVTKPoints::setUp();
                _data = new TestDataWriterVTKPoints_Data();CPPUNIT_ASSERT(_data);

                _data->timestepFilename = "hex8_points.vtk";
                _data->vertexFilename = "hex8_points_vertex.vtk";

                TestDataWriterPoints::_setDataHex();
                TestDataWriterPoints::_initialize();

                PYLITH_METHOD_END;
            } // setUp

        }; // class TestDataWriterVTKPoints_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDataWriterVTKPoints_Hex);

    } // meshio
} // pylith

// End of file
