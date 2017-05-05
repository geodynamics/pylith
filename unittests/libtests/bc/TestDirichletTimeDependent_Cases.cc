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

#include "TestDirichletTimeDependent.hh" // Implementation of cases

// ----------------------------------------------------------------------
namespace pylith {
    namespace bc {

        // --------------------------------------------------------------
        class TestDirichletTimeDependent_TriP1 : public TestDirichletTimeDependent {
            CPPUNIT_TEST_SUB_SUITE(TestDirichletTimeDependent_TriP1, TestDirichletTimeDependent);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestDirichletTimeDependent::setUp();
                _data = new TestDirichletTimeDependent_Data();CPPUNIT_ASSERT(_data);
                _data->meshFilename = "data/tri_small.mesh";
                _data->bcLabel = "boundary_bottom";

                _data->lengthScale = 1000.0;
                _data->timeScale = 10.0;
                _data->pressureScale = 0.1;
                _data->densityScale = 1.0;

                _data->field = "displacement";
                _data->isConstrainedFieldVector = true;
                _data->numConstrainedDOF = 1;
                static const int constrainedDOF[1] = { 1 };
                _data->constrainedDOF = const_cast<int*>(constrainedDOF);

                _data->useInitial = true;
                _data->useRate = false;
                _data->useTimeHistory = false;

                _data->numSolnFields = 2;
                static const char* solutionFields[2] = { "displacement", "velocity" };
                _data->solutionFields = const_cast<const char**>(solutionFields);
                static const pylith::topology::Field::DiscretizeInfo solnDiscretizations[2] = {
                    {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // displacement
                    {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // velocity
                };
                _data->solnDiscretizations = const_cast<pylith::topology::Field::DiscretizeInfo*>(solnDiscretizations);
                // solnDiscretizations

                _data->numAuxFields = 1;
                static const char* auxFields[1] = { "initial_amplitude" };
                _data->auxFields = const_cast<const char**>(auxFields);
                static const pylith::topology::Field::DiscretizeInfo auxDiscretizations[1] = {
                    {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // initial_amplitude
                };
                _data->auxDiscretizations = const_cast<pylith::topology::Field::DiscretizeInfo*>(auxDiscretizations);
                // auxDBFilename

                // t
                // solnDBFilename
            } // setUp
        }; // class TestDirichletTimeDependent_TriP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDirichletTimeDependent_TriP1);

#if 0
        // --------------------------------------------------------------
        class TestDirichletTimeDependent_QuadP1 : public TestDirichletTimeDependent {
            CPPUNIT_TEST_SUB_SUITE(TestDirichletTimeDependent_QuadP1, TestDirichletTimeDependent);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestDirichletTimeDependent::setUp();
                _data = new TestDirichletTimeDependent_Data();CPPUNIT_ASSERT(_data);
                _data->meshFilename = "data/quad_small.mesh";
                _data->bcLabel = "boundary_top";
            } // setUp
        }; // class TestDirichletTimeDependent_QuadP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDirichletTimeDependent_QuadP1);


        // --------------------------------------------------------------
        class TestDirichletTimeDependent_TetP1 : public TestDirichletTimeDependent {
            CPPUNIT_TEST_SUB_SUITE(TestDirichletTimeDependent_TetP1, TestDirichletTimeDependent);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestDirichletTimeDependent::setUp();
                _data = new TestDirichletTimeDependent_Data();CPPUNIT_ASSERT(_data);
            } // setUp
        }; // class TestDirichletTimeDependent_TetP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDirichletTimeDependent_TetP1);


        // --------------------------------------------------------------
        class TestDirichletTimeDependent_Hex : public TestDirichletTimeDependent {
            CPPUNIT_TEST_SUB_SUITE(TestDirichletTimeDependent_Hex, TestDirichletTimeDependent);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestDirichletTimeDependent::setUp();
                _data = new TestDirichletTimeDependent_Data();CPPUNIT_ASSERT(_data);
            } // setUp
        }; // class TestDirichletTimeDependent_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDirichletTimeDependent_Hex);
#endif

    } // namespace bc
} // namespace pylith


// End of file
