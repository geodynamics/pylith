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
// Copyright (c) 2010-2022 University of California, Davis
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

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace problems {
        class TestSolutionFactory_Cases;
    } // problems
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::problems::TestSolutionFactory_Cases {
public:

    // Data factory methods
    static TestSolutionFactory_Data* Tri(void);

    static TestSolutionFactory_Data* Hex(void);

private:

    static
    double displacement_2d_x(const double x,
                             const double y) {
        return 3.0*x + 2.0*y;
    } // displacement_x

    static
    double displacement_2d_y(const double x,
                             const double y) {
        return 2.0*x - 0.1*y;
    } // displacement_y

    static
    double displacement_3d_x(const double x,
                             const double y,
                             const double z) {
        return 3.0*x + 2.0*y - 0.4*z;
    } // displacement_x

    static
    double displacement_3d_y(const double x,
                             const double y,
                             const double z) {
        return 2.0*x - 0.1*y + 0.2*z;
    } // displacement_y

    static
    double displacement_3d_z(const double x,
                             const double y,
                             const double z) {
        return 0.3*x + 0.8*y + 0.1*z;
    } // displacement_z

    static
    const char* displacement_units(void) {
        return "m";
    } // displacement_units

    static
    double pressure_2d(const double x,
                       const double y) {
        return -0.3*x*y + 0.2*y*y;
    } // pressure

    static
    double pressure_3d(const double x,
                       const double y,
                       const double z) {
        return -0.3*x*y + 0.2*y*y + 0.7*y*z;
    } // pressure

    static
    const char* pressure_units(void) {
        return "Pa";
    } // pressure_units

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestSolutionFactory::Tri::testDispVel", "[TestSolutionFactory][testDispVel]") {
    pylith::problems::TestSolutionFactory(pylith::problems::TestSolutionFactory_Cases::Tri()).testDispVel();
}
TEST_CASE("TestSolutionFactory::Tri::testDispLagrangeFault", "[TestSolutionFactory][testDispLagrangeFault]") {
    pylith::problems::TestSolutionFactory(pylith::problems::TestSolutionFactory_Cases::Tri()).testDispLagrangeFault();
}
TEST_CASE("TestSolutionFactory::Tri::testPressure", "[TestSolutionFactory][testPressure]") {
    pylith::problems::TestSolutionFactory(pylith::problems::TestSolutionFactory_Cases::Tri()).testPressure();
}
TEST_CASE("TestSolutionFactory::Tri::testDispTemp", "[TestSolutionFactory][testDispTemp]") {
    pylith::problems::TestSolutionFactory(pylith::problems::TestSolutionFactory_Cases::Tri()).testDispTemp();
}
TEST_CASE("TestSolutionFactory::Tri::testSetValues", "[TestSolutionFactory][testSetValues]") {
    pylith::problems::TestSolutionFactory(pylith::problems::TestSolutionFactory_Cases::Tri()).testSetValues();
}

TEST_CASE("TestSolutionFactory::Hex::testDispVel", "[TestSolutionFactory][testDispVel]") {
    pylith::problems::TestSolutionFactory(pylith::problems::TestSolutionFactory_Cases::Hex()).testDispVel();
}
TEST_CASE("TestSolutionFactory::Hex::testDispLagrangeFault", "[TestSolutionFactory][testDispLagrangeFault]") {
    pylith::problems::TestSolutionFactory(pylith::problems::TestSolutionFactory_Cases::Hex()).testDispLagrangeFault();
}
TEST_CASE("TestSolutionFactory::Hex::testPressure", "[TestSolutionFactory][testPressure]") {
    pylith::problems::TestSolutionFactory(pylith::problems::TestSolutionFactory_Cases::Hex()).testPressure();
}
TEST_CASE("TestSolutionFactory::Hex::testDispTemp", "[TestSolutionFactory][testDispTemp]") {
    pylith::problems::TestSolutionFactory(pylith::problems::TestSolutionFactory_Cases::Hex()).testDispTemp();
}
TEST_CASE("TestSolutionFactory::Hex::testSetValues", "[TestSolutionFactory][testSetValues]") {
    pylith::problems::TestSolutionFactory(pylith::problems::TestSolutionFactory_Cases::Hex()).testSetValues();
}

// ------------------------------------------------------------------------------------------------
pylith::problems::TestSolutionFactory_Data*
pylith::problems::TestSolutionFactory_Cases::Tri(void) {
    pylith::problems::TestSolutionFactory_Data* data = new pylith::problems::TestSolutionFactory_Data();
    assert(data);

    data->dimension = 2;
    data->meshFilename = "data/tri.mesh";
    data->cs = new spatialdata::geocoords::CSCart();assert(data->cs);
    data->cs->setSpaceDim(data->dimension);

    assert(data->solutionDB);
    data->solutionDB->addValue("displacement_x", displacement_2d_x, displacement_units());
    data->solutionDB->addValue("displacement_y", displacement_2d_y, displacement_units());
    data->solutionDB->addValue("pressure", pressure_2d, pressure_units());
    data->solutionDB->setDescription("solution");
    data->solutionDB->setCoordSys(*data->cs);

    data->subfields["displacement"].description.numComponents = 2;
    data->subfields["velocity"].description.numComponents = 2;
    data->subfields["lagrange_multiplier_fault"].description.numComponents = 2;

    return data;
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::problems::TestSolutionFactory_Data*
pylith::problems::TestSolutionFactory_Cases::Hex(void) {
    pylith::problems::TestSolutionFactory_Data* data = new pylith::problems::TestSolutionFactory_Data();
    assert(data);

    data->dimension = 3;
    data->meshFilename = "data/hex.mesh";
    data->cs = new spatialdata::geocoords::CSCart();assert(data->cs);
    data->cs->setSpaceDim(data->dimension);

    assert(data->solutionDB);
    data->solutionDB->addValue("displacement_x", displacement_3d_x, displacement_units());
    data->solutionDB->addValue("displacement_y", displacement_3d_y, displacement_units());
    data->solutionDB->addValue("displacement_z", displacement_3d_z, displacement_units());
    data->solutionDB->addValue("pressure", pressure_3d, pressure_units());
    data->solutionDB->setDescription("solution");
    data->solutionDB->setCoordSys(*data->cs);

    return data;
} // Hex


// End of file
