// -*- C++ -*-
//
// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
//

#include <portinfo>

#include "TestReverseCuthillMcKee.hh" // Implementation of class methods

// -----------------------------------------------------------------------------
namespace pylith {
    namespace topology {

        // ---------------------------------------------------------------------
        class TestReverseCuthillMcKee_Tri_Nofault : public TestReverseCuthillMcKee {

            CPPUNIT_TEST_SUB_SUITE( TestReverseCuthillMcKee_Tri_Nofault, TestReverseCuthillMcKee );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestReverseCuthillMcKee::setUp();

                _data->filename = "data/reorder_tri3.mesh";
                _data->faultLabel = NULL;
            }   // setUp


        };  // TestReverseCuthillMcKee_Tri_Nofault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestReverseCuthillMcKee_Tri_Nofault );

        // ---------------------------------------------------------------------
        class TestReverseCuthillMcKee_Tri_Fault : public TestReverseCuthillMcKee {

            CPPUNIT_TEST_SUB_SUITE( TestReverseCuthillMcKee_Tri_Fault, TestReverseCuthillMcKee );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestReverseCuthillMcKee::setUp();

                _data->filename = "data/reorder_tri3.mesh";
                _data->faultLabel = "fault";
            }   // setUp


        };  // TestReverseCuthillMcKee_Tri_Fault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestReverseCuthillMcKee_Tri_Fault );

        // ---------------------------------------------------------------------
        class TestReverseCuthillMcKee_Quad_Nofault : public TestReverseCuthillMcKee {

            CPPUNIT_TEST_SUB_SUITE( TestReverseCuthillMcKee_Quad_Nofault, TestReverseCuthillMcKee );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestReverseCuthillMcKee::setUp();

                _data->filename = "data/reorder_quad4.mesh";
                _data->faultLabel = NULL;
            }   // setUp


        };  // TestReverseCuthillMcKee_Quad_Nofault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestReverseCuthillMcKee_Quad_Nofault );

        // ---------------------------------------------------------------------
        class TestReverseCuthillMcKee_Quad_Fault : public TestReverseCuthillMcKee {

            CPPUNIT_TEST_SUB_SUITE( TestReverseCuthillMcKee_Quad_Fault, TestReverseCuthillMcKee );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestReverseCuthillMcKee::setUp();

                _data->filename = "data/reorder_quad4.mesh";
                _data->faultLabel = "fault";
            }   // setUp


        };  // TestReverseCuthillMcKee_Quad_Fault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestReverseCuthillMcKee_Quad_Fault );

        // ---------------------------------------------------------------------
        class TestReverseCuthillMcKee_Tet_Nofault : public TestReverseCuthillMcKee {

            CPPUNIT_TEST_SUB_SUITE( TestReverseCuthillMcKee_Tet_Nofault, TestReverseCuthillMcKee );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestReverseCuthillMcKee::setUp();

                _data->filename = "data/reorder_tet4.mesh";
                _data->faultLabel = NULL;
            }   // setUp


        };  // TestReverseCuthillMcKee_Tet_Nofault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestReverseCuthillMcKee_Tet_Nofault );

        // ---------------------------------------------------------------------
        class TestReverseCuthillMcKee_Tet_Fault : public TestReverseCuthillMcKee {

            CPPUNIT_TEST_SUB_SUITE( TestReverseCuthillMcKee_Tet_Fault, TestReverseCuthillMcKee );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestReverseCuthillMcKee::setUp();

                _data->filename = "data/reorder_tet4.mesh";
                _data->faultLabel = "fault";
            }   // setUp


        };  // TestReverseCuthillMcKee_Tet_Fault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestReverseCuthillMcKee_Tet_Fault );

        // ---------------------------------------------------------------------
        class TestReverseCuthillMcKee_Hex_Nofault : public TestReverseCuthillMcKee {

            CPPUNIT_TEST_SUB_SUITE( TestReverseCuthillMcKee_Hex_Nofault, TestReverseCuthillMcKee );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestReverseCuthillMcKee::setUp();

                _data->filename = "data/reorder_hex8.mesh";
                _data->faultLabel = NULL;
            }   // setUp


        };  // TestReverseCuthillMcKee_Hex_Nofault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestReverseCuthillMcKee_Hex_Nofault );

        // ---------------------------------------------------------------------
        class TestReverseCuthillMcKee_Hex_Fault : public TestReverseCuthillMcKee {

            CPPUNIT_TEST_SUB_SUITE( TestReverseCuthillMcKee_Hex_Fault, TestReverseCuthillMcKee );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestReverseCuthillMcKee::setUp();

                _data->filename = "data/reorder_hex8.mesh";
                _data->faultLabel = "fault";
            }   // setUp


        };  // TestReverseCuthillMcKee_Hex_Fault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestReverseCuthillMcKee_Hex_Fault );

    }   // topology
}   // pylith


// End of file
