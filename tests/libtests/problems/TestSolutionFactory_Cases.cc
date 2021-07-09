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

#include "TestSolutionFactory.hh" // Implementation of cases

#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

namespace pylith {
    namespace problems {
        class TestSolutionFactory_Tri;
        class TestSolutionFactory_Hex;
    } // problems
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::problems::TestSolutionFactory_Tri : public pylith::problems::TestSolutionFactory {
    CPPUNIT_TEST_SUB_SUITE(TestSolutionFactory_Tri, TestSolutionFactory);
    CPPUNIT_TEST_SUITE_END();

    static
    double displacement_x(const double x,
                          const double y) {
        return 3.0*x + 2.0*y;
    } // displacement_x

    static
    double displacement_y(const double x,
                          const double y) {
        return 2.0*x - 0.1*y;
    } // displacement_y

    static
    const char* displacement_units(void) {
        return "m";
    } // displacement_units

    static
    double pressure(const double x,
                    const double y) {
        return -0.3*x*y + 0.2*y*y;
    } // pressure

    static
    const char* pressure_units(void) {
        return "Pa";
    } // pressure_units

protected:

    void setUp(void) {
        TestSolutionFactory::setUp();

        CPPUNIT_ASSERT(_data);
        _data->dimension = 2;
        _data->meshFilename = "data/tri.mesh";
        _data->cs = new spatialdata::geocoords::CSCart();CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->dimension);

        CPPUNIT_ASSERT(_data->solutionDB);
        _data->solutionDB->addValue("displacement_x", displacement_x, displacement_units());
        _data->solutionDB->addValue("displacement_y", displacement_y, displacement_units());
        _data->solutionDB->addValue("pressure", pressure, pressure_units());
        _data->solutionDB->setLabel("solution");
        _data->solutionDB->setCoordSys(*_data->cs);

        _data->subfields["displacement"].description.numComponents = 2;
        _data->subfields["velocity"].description.numComponents = 2;
        _data->subfields["lagrange_multiplier_fault"].description.numComponents = 2;

        _initialize();
    } // setUp

}; // TestSolutionFactory_Tri
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::problems::TestSolutionFactory_Tri);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::problems::TestSolutionFactory_Hex : public pylith::problems::TestSolutionFactory {
    CPPUNIT_TEST_SUB_SUITE(TestSolutionFactory_Hex, TestSolutionFactory);
    CPPUNIT_TEST_SUITE_END();

    static
    double displacement_x(const double x,
                          const double y,
                          const double z) {
        return 3.0*x + 2.0*y - 0.4*z;
    } // displacement_x

    static
    double displacement_y(const double x,
                          const double y,
                          const double z) {
        return 2.0*x - 0.1*y + 0.2*z;
    } // displacement_y

    static
    double displacement_z(const double x,
                          const double y,
                          const double z) {
        return 0.3*x + 0.8*y + 0.1*z;
    } // displacement_z

    static
    const char* displacement_units(void) {
        return "m";
    } // displacement_units

    static
    double pressure(const double x,
                    const double y,
                    const double z) {
        return -0.3*x*y + 0.2*y*y + 0.7*y*z;
    } // pressure

    static
    const char* pressure_units(void) {
        return "Pa";
    } // pressure_units

protected:

    void setUp(void) {
        TestSolutionFactory::setUp();

        CPPUNIT_ASSERT(_data);
        _data->dimension = 3;
        _data->meshFilename = "data/hex.mesh";
        _data->cs = new spatialdata::geocoords::CSCart();CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->dimension);

        CPPUNIT_ASSERT(_data->solutionDB);
        _data->solutionDB->addValue("displacement_x", displacement_x, displacement_units());
        _data->solutionDB->addValue("displacement_y", displacement_y, displacement_units());
        _data->solutionDB->addValue("displacement_z", displacement_z, displacement_units());
        _data->solutionDB->addValue("pressure", pressure, pressure_units());
        _data->solutionDB->setLabel("solution");
        _data->solutionDB->setCoordSys(*_data->cs);

        _initialize();
    } // setUp

}; // TestSolutionFactory_Hex
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::problems::TestSolutionFactory_Hex);

// End of file
