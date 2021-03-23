// -*- C++ -*-
//
// -----------------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// -----------------------------------------------------------------------------
//

#include <portinfo>

#include "TestFaultCohesive.hh" // Implementation of class methods

// -----------------------------------------------------------------------------
namespace pylith {
    namespace faults {
        // ---------------------------------------------------------------------
        class TestFaultCohesive_Quad : public TestFaultCohesive {
            CPPUNIT_TEST_SUB_SUITE(TestFaultCohesive_Quad, TestFaultCohesive);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestFaultCohesive::setUp();

                _data->filename = "data/quad_patches_a.mesh";
                _data->faultLabel = "fault";
                _data->edgeLabel = NULL;
            } // setUp

        }; // TestFaultCohesive_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestFaultCohesive_Quad);

    } // faults
} // pylith

// End of file
