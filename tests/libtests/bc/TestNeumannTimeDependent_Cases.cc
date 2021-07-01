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

#include "TestNeumannTimeDependent.hh" // Implementation of cases

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
namespace pylith {
    namespace bc {
        // --------------------------------------------------------------
        class TestNeumannTimeDependent_TriP1 : public TestNeumannTimeDependent {
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
            static double initial_traction_tangential(const double x,
                                                      const double y) {
                return -1.5*x + 2.4*y;
            } // initial_traction_tangential

            static double initial_traction_normal(const double x,
                                                  const double y) {
                return 2.4*x + 1.8*y;
            } // initial_traction_normal

            // Solution field at time t.

            static double disp_x(const double x,
                                 const double y) {
                return FILL_VALUE;
            } // disp_x

            static double disp_y(const double x,
                                 const double y) {
                return FILL_VALUE;
            } // disp_y

            static double vel_x(const double x,
                                const double y) {
                return FILL_VALUE;
            } // vel_x

            static double vel_y(const double x,
                                const double y) {
                return FILL_VALUE;
            } // vel_y

            static double fluid_press(const double x,
                                      const double y) {
                return FILL_VALUE;
            } // fluid_press

            CPPUNIT_TEST_SUB_SUITE(TestNeumannTimeDependent_TriP1, TestNeumannTimeDependent);
            CPPUNIT_TEST_SUITE_END();

protected:

            void setUp(void) {
                TestNeumannTimeDependent::setUp();
                _data = new TestNeumannTimeDependent_Data();CPPUNIT_ASSERT(_data);
                _data->meshFilename = "data/tri_small.mesh";
                _data->bcLabel = "boundary_bottom";
                _data->cs = new spatialdata::geocoords::CSCart();CPPUNIT_ASSERT(_data->cs);
                _data->cs->setSpaceDim(2);

                CPPUNIT_ASSERT(_data->normalizer);
                _data->normalizer->setLengthScale(1000.0);
                _data->normalizer->setTimeScale(10.0);
                _data->normalizer->setPressureScale(0.1);
                _data->normalizer->setDensityScale(2.0);

                _data->field = "displacement";
                _data->vectorFieldType = pylith::topology::Field::VECTOR;
                _data->scale = _data->normalizer->getPressureScale();

                _data->useInitial = true;
                _data->useRate = false;
                _data->useTimeHistory = false;

                _data->numAuxSubfields = 1;
                static const char* auxSubfields[1] = { "initial_amplitude" };
                _data->auxSubfields = const_cast<const char**>(auxSubfields);
                static const pylith::topology::Field::Discretization auxDiscretizations[1] = {
                    pylith::topology::Field::Discretization(1, 1), // initial_amplitude
                };
                _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(auxDiscretizations);

                CPPUNIT_ASSERT(_data->auxDB);
                _data->auxDB->setCoordSys(*_data->cs);
                _data->auxDB->addValue("initial_amplitude_tangential", initial_traction_tangential, pressure_units());
                _data->auxDB->addValue("initial_amplitude_normal", initial_traction_normal, pressure_units());

                _data->t = 1.23;
                _data->solnNumSubfields = 3;
                static const pylith::topology::Field::Discretization solnDiscretizations[3] = {
                    pylith::topology::Field::Discretization(1, 1), // displacement
                    pylith::topology::Field::Discretization(1, 1), // velocity
                    pylith::topology::Field::Discretization(1, 1), // fluid_pressure
                };
                _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(solnDiscretizations);

                CPPUNIT_ASSERT(_data->solnDB);
                _data->solnDB->setCoordSys(*_data->cs);
                _data->solnDB->addValue("displacement_x", disp_x, disp_units());
                _data->solnDB->addValue("displacement_y", disp_y, disp_units());
                _data->solnDB->addValue("velocity_x", vel_x, vel_units());
                _data->solnDB->addValue("velocity_y", vel_y, vel_units());
                _data->solnDB->addValue("fluid_pressure", fluid_press, pressure_units());

            } // setUp

        }; // class TestNeumannTimeDependent_TriP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestNeumannTimeDependent_TriP1);

#if 0
        // --------------------------------------------------------------
        class TestNeumannTimeDependent_QuadP1 : public TestNeumannTimeDependent {
            CPPUNIT_TEST_SUB_SUITE(TestNeumannTimeDependent_QuadP1, TestNeumannTimeDependent);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestNeumannTimeDependent::setUp();
                _data = new TestNeumannTimeDependent_Data();CPPUNIT_ASSERT(_data);
                _data->meshFilename = "data/quad_small.mesh";
                _data->bcLabel = "boundary_top";
            } // setUp

        }; // class TestNeumannTimeDependent_QuadP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestNeumannTimeDependent_QuadP1);

        // --------------------------------------------------------------
        class TestNeumannTimeDependent_TetP1 : public TestNeumannTimeDependent {
            CPPUNIT_TEST_SUB_SUITE(TestNeumannTimeDependent_TetP1, TestNeumannTimeDependent);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestNeumannTimeDependent::setUp();
                _data = new TestNeumannTimeDependent_Data();CPPUNIT_ASSERT(_data);
            } // setUp

        }; // class TestNeumannTimeDependent_TetP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestNeumannTimeDependent_TetP1);

        // --------------------------------------------------------------
        class TestNeumannTimeDependent_Hex : public TestNeumannTimeDependent {
            CPPUNIT_TEST_SUB_SUITE(TestNeumannTimeDependent_Hex, TestNeumannTimeDependent);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestNeumannTimeDependent::setUp();
                _data = new TestNeumannTimeDependent_Data();CPPUNIT_ASSERT(_data);
            } // setUp

        }; // class TestNeumannTimeDependent_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestNeumannTimeDependent_Hex);
#endif

    } // namespace bc
} // namespace pylith

// End of file
