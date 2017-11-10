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

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional


// ----------------------------------------------------------------------
namespace pylith {
    namespace bc {

        // --------------------------------------------------------------
        class TestDirichletTimeDependent_TriP1 : public TestDirichletTimeDependent {

            static const char* disp_units(void) {
                return "m";
            }
            static const char* vel_units(void) {
                return "m/s";
            }
            static const char* pressure_units(void) {
                return "Pa";
            }

            /// Spatial database user functions for auxiliary subfields.

            // Initial amplitude
            static double initial_disp_x(const double x,
                                         const double y) {
                return -888.0;
            } // initial_disp_x
            static double initial_disp_y(const double x,
                                         const double y) {
                return 2.4*x + 1.8*y;
            } // initial_disp_y

            // Solution field at time t.

            static double disp_x(const double x,
                                 const double y) {
                return -999.0;
            } // disp_x
            static double disp_y(const double x,
                                 const double y) {
                return initial_disp_y(x,y);
            } // disp_y
            static double vel_x(const double x,
                                const double y) {
                return -999.0;
            } // vel_x
            static double vel_y(const double x,
                                const double y) {
                return -999.0;
            } // vel_y
            static double fluid_press(const double x,
                                      const double y) {
                return -999.0;
            } // fluid_press

            CPPUNIT_TEST_SUB_SUITE(TestDirichletTimeDependent_TriP1, TestDirichletTimeDependent);
            CPPUNIT_TEST_SUITE_END();

protected:
            void setUp(void) {
                TestDirichletTimeDependent::setUp();
                _data = new TestDirichletTimeDependent_Data();CPPUNIT_ASSERT(_data);
                _data->meshFilename = "data/tri_small.mesh";
                _data->bcLabel = "boundary_bottom";
                _data->cs = new spatialdata::geocoords::CSCart();CPPUNIT_ASSERT(_data->cs);
                _data->cs->setSpaceDim(2);
                _data->cs->initialize();

                CPPUNIT_ASSERT(_data->normalizer);
                _data->normalizer->lengthScale(1000.0);
                _data->normalizer->timeScale(10.0);
                _data->normalizer->pressureScale(0.1);
                _data->normalizer->densityScale(2.0);

                _data->field = "displacement";
                _data->vectorFieldType = pylith::topology::Field::VECTOR;
                _data->numConstrainedDOF = 1;
                static const int constrainedDOF[1] = { 1 };
                _data->constrainedDOF = const_cast<int*>(constrainedDOF);

                _data->useInitial = true;
                _data->useRate = false;
                _data->useTimeHistory = false;

                _data->numAuxSubfields = 1;
                static const char* auxSubfields[1] = { "initial_amplitude" };
                _data->auxSubfields = const_cast<const char**>(auxSubfields);
                static const pylith::topology::Field::Discretization auxDiscretizations[1] = {
                    {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // initial_amplitude
                };
                _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(auxDiscretizations);

                CPPUNIT_ASSERT(_data->auxDB);
                _data->auxDB->coordsys(*_data->cs);
                _data->auxDB->addValue("initial_amplitude_x", initial_disp_x, disp_units());
                _data->auxDB->addValue("initial_amplitude_y", initial_disp_y, disp_units());

                _data->t = 1.23;
                _data->solnNumSubfields = 3;
                static const pylith::topology::Field::Discretization solnDiscretizations[3] = {
                    {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // displacement
                    {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // velocity
                    {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // fluid_pressure
                };
                _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(solnDiscretizations);

                CPPUNIT_ASSERT(_data->solnDB);
                _data->solnDB->coordsys(*_data->cs);
                _data->solnDB->addValue("displacement_x", disp_x, disp_units());
                _data->solnDB->addValue("displacement_y", disp_y, disp_units());
                _data->solnDB->addValue("velocity_x", vel_x, vel_units());
                _data->solnDB->addValue("velocity_y", vel_y, vel_units());
                _data->solnDB->addValue("fluid_pressure", fluid_press, pressure_units());

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
