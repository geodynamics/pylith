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

#include "TestXdmf.hh" // Implementation of class methods

// ----------------------------------------------------------------------
namespace pylith {
    namespace meshio {

        // --------------------------------------------------------------
        class TestXdmf_TriP1_Vertex : public TestXdmf {
            CPPUNIT_TEST_SUB_SUITE(TestXdmf_TriP1_Vertex, TestXdmf);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestXdmf::setUp();
                _data->filenameHDF5 = "data/tri3_vertex.h5";
                _data->filenameXdmf = "tri3_vertex.xmf";
            } // setUp
        }; // class TestXdmf_TriP1_Vertex
        CPPUNIT_TEST_SUITE_REGISTRATION( TestXdmf_TriP1_Vertex );


        // --------------------------------------------------------------
        class TestXdmf_TriP1_Cell : public TestXdmf {
            CPPUNIT_TEST_SUB_SUITE(TestXdmf_TriP1_Cell, TestXdmf);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestXdmf::setUp();
                _data->filenameHDF5 = "data/tri3_cell.h5";
                _data->filenameXdmf = "tri3_cell.xmf";
            } // setUp
        }; // class TestXdmf_TriP1_Cell
        CPPUNIT_TEST_SUITE_REGISTRATION( TestXdmf_TriP1_Cell );

        // --------------------------------------------------------------
        class TestXdmf_QuadP1_Vertex : public TestXdmf {
            CPPUNIT_TEST_SUB_SUITE(TestXdmf_QuadP1_Vertex, TestXdmf);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestXdmf::setUp();
                _data->filenameHDF5 = "data/quad4_vertex.h5";
                _data->filenameXdmf = "quad4_vertex.xmf";
            } // setUp
        }; // class TestXdmf_QuadP1_Vertex
        CPPUNIT_TEST_SUITE_REGISTRATION( TestXdmf_QuadP1_Vertex );


        // --------------------------------------------------------------
        class TestXdmf_QuadP1_Cell : public TestXdmf {
            CPPUNIT_TEST_SUB_SUITE(TestXdmf_QuadP1_Cell, TestXdmf);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestXdmf::setUp();
                _data->filenameHDF5 = "data/quad4_cell.h5";
                _data->filenameXdmf = "quad4_cell.xmf";
            } // setUp
        }; // class TestXdmf_QuadP1_Cell
        CPPUNIT_TEST_SUITE_REGISTRATION( TestXdmf_QuadP1_Cell );

        // --------------------------------------------------------------
        class TestXdmf_TetP1_Vertex : public TestXdmf {
            CPPUNIT_TEST_SUB_SUITE(TestXdmf_TetP1_Vertex, TestXdmf);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestXdmf::setUp();
                _data->filenameHDF5 = "data/tet4_vertex.h5";
                _data->filenameXdmf = "tet4_vertex.xmf";
            } // setUp
        }; // class TestXdmf_TetP1_Vertex
        CPPUNIT_TEST_SUITE_REGISTRATION( TestXdmf_TetP1_Vertex );


        // --------------------------------------------------------------
        class TestXdmf_TetP1_Cell : public TestXdmf {
            CPPUNIT_TEST_SUB_SUITE(TestXdmf_TetP1_Cell, TestXdmf);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestXdmf::setUp();
                _data->filenameHDF5 = "data/tet4_cell.h5";
                _data->filenameXdmf = "tet4_cell.xmf";
            } // setUp
        }; // class TestXdmf_TetP1_Cell
        CPPUNIT_TEST_SUITE_REGISTRATION( TestXdmf_TetP1_Cell );

        // --------------------------------------------------------------
        class TestXdmf_HexP1_Vertex : public TestXdmf {
            CPPUNIT_TEST_SUB_SUITE(TestXdmf_HexP1_Vertex, TestXdmf);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestXdmf::setUp();
                _data->filenameHDF5 = "data/hex8_vertex.h5";
                _data->filenameXdmf = "hex8_vertex.xmf";
            } // setUp
        }; // class TestXdmf_HexP1_Vertex
        CPPUNIT_TEST_SUITE_REGISTRATION( TestXdmf_HexP1_Vertex );


        // --------------------------------------------------------------
        class TestXdmf_HexP1_Cell : public TestXdmf {
            CPPUNIT_TEST_SUB_SUITE(TestXdmf_HexP1_Cell, TestXdmf);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestXdmf::setUp();
                _data->filenameHDF5 = "data/hex8_cell.h5";
                _data->filenameXdmf = "hex8_cell.xmf";
            } // setUp
        }; // class TestXdmf_HexP1_Cell
        CPPUNIT_TEST_SUITE_REGISTRATION( TestXdmf_HexP1_Cell );

    } // meshio
} // pylith


// End of file 
