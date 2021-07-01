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

#include "TestAbsorbingDampers.hh" // Implementation of cases

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
namespace pylith {
    namespace bc {
        // --------------------------------------------------------------
        class TestAbsorbingDampers_TriP1 : public TestAbsorbingDampers {
            static const char* density_units(void) {
                return "kg/m**3";
            } // density_units

            static const char* disp_units(void) {
                return "m";
            } // disp_units

            static const char* vel_units(void) {
                return "m/s";
            } // vel_units

            static const char* pressure_units(void) {
                return "Pa";
            } // pressure_units

            /// Spatial database user functions for auxiliary subfields.

            static double density(const double x,
                                  const double y) {
                return 2500.0;
            } // density

            static double vs(const double x,
                             const double y) {
                return 3000.0;
            } // vs

            static double vp(const double x,
                             const double y) {
                return vs(x,y)*sqrt(3.0);
            } // vp

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
                return FILL_VALUE; // :TODO:
            } // vel_x

            static double vel_y(const double x,
                                const double y) {
                return FILL_VALUE; // :TODO:
            } // vel_y

            static double fluid_press(const double x,
                                      const double y) {
                return FILL_VALUE;
            } // fluid_press

            CPPUNIT_TEST_SUB_SUITE(TestAbsorbingDampers_TriP1, TestAbsorbingDampers);
            CPPUNIT_TEST_SUITE_END();

protected:

            void setUp(void) {
                TestAbsorbingDampers::setUp();
                _data = new TestAbsorbingDampers_Data();CPPUNIT_ASSERT(_data);
                _data->meshFilename = "data/tri_small.mesh";
                _data->bcLabel = "boundary_bottom";
                _data->cs = new spatialdata::geocoords::CSCart();CPPUNIT_ASSERT(_data->cs);
                _data->cs->setSpaceDim(2);

                CPPUNIT_ASSERT(_data->normalizer);
                _data->normalizer->setLengthScale(1000.0);
                _data->normalizer->setTimeScale(10.0);
                _data->normalizer->setPressureScale(0.1);
                _data->normalizer->setDensityScale(2.0);

                _data->field = "velocity";
                _data->vectorFieldType = pylith::topology::Field::VECTOR;

                _data->numAuxSubfields = 3;
                static const char* auxSubfields[3] = {
                    "density",
                    "vp",
                    "vs",
                };
                _data->auxSubfields = const_cast<const char**>(auxSubfields);
                static const pylith::topology::Field::Discretization auxDiscretizations[3] = {
                    pylith::topology::Field::Discretization(0, 1), // density
                    pylith::topology::Field::Discretization(0, 1), // vp
                    pylith::topology::Field::Discretization(0, 1), // vs
                };
                _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(auxDiscretizations);

                CPPUNIT_ASSERT(_data->auxDB);
                _data->auxDB->setCoordSys(*_data->cs);
                _data->auxDB->addValue("density", density, density_units());
                _data->auxDB->addValue("vp", vp, vel_units());
                _data->auxDB->addValue("vs", vs, vel_units());

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

        }; // class TestAbsorbingDampers_TriP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAbsorbingDampers_TriP1);

#if 0
        // --------------------------------------------------------------
        class TestAbsorbingDampers_QuadP1 : public TestAbsorbingDampers {
            CPPUNIT_TEST_SUB_SUITE(TestAbsorbingDampers_QuadP1, TestAbsorbingDampers);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAbsorbingDampers::setUp();
                _data = new TestAbsorbingDampers_Data();CPPUNIT_ASSERT(_data);
                _data->meshFilename = "data/quad_small.mesh";
                _data->bcLabel = "boundary_top";
            } // setUp

        }; // class TestAbsorbingDampers_QuadP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAbsorbingDampers_QuadP1);

        // --------------------------------------------------------------
        class TestAbsorbingDampers_TetP1 : public TestAbsorbingDampers {
            CPPUNIT_TEST_SUB_SUITE(TestAbsorbingDampers_TetP1, TestAbsorbingDampers);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAbsorbingDampers::setUp();
                _data = new TestAbsorbingDampers_Data();CPPUNIT_ASSERT(_data);
            } // setUp

        }; // class TestAbsorbingDampers_TetP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAbsorbingDampers_TetP1);

        // --------------------------------------------------------------
        class TestAbsorbingDampers_Hex : public TestAbsorbingDampers {
            CPPUNIT_TEST_SUB_SUITE(TestAbsorbingDampers_Hex, TestAbsorbingDampers);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAbsorbingDampers::setUp();
                _data = new TestAbsorbingDampers_Data();CPPUNIT_ASSERT(_data);
            } // setUp

        }; // class TestAbsorbingDampers_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAbsorbingDampers_Hex);
#endif

    } // namespace bc
} // namespace pylith

// End of file
