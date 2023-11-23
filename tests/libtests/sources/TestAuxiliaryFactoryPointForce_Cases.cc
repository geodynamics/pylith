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

#include "TestAuxiliaryFactoryPointForce.hh" // Implementation of cases

#include "pylith/sources/AuxiliaryFactoryPointForce.hh" // USES AuxiliaryFactoryPointForce
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "catch2/catch_test_macros.hpp"

#include <cmath> // USES fabs()

// forward declarations
namespace pylith {
    namespace sources {
        class TestAuxiliaryFactoryPointForce_Cases;
    } // sources
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::sources::TestAuxiliaryFactoryPointForce_Cases {
public:

    // Data factory methods
    static TestAuxiliaryFactoryPointForce_Data* Tri(void);

    static TestAuxiliaryFactoryPointForce_Data* Hex(void);

private:

    static
    double point_force_2d_x(const double x,
                            const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // point_force_2d_x

    static
    double point_force_2d_y(const double x,
                            const double y) {
        return -0.8*x*y + 0.2*y*y;
    } // point_force_2d_yy

    static
    double point_force_2d_z(const double x,
                            const double y) {
        return -0.3*x*x + 0.5*y*y;
    } // point_force_2d_z

    static
    const char* point_force_units(void) {
        return "Pa";
    } // point_force_units

    static
    double time_delay_2d(const double x,
                         const double y) {
        return 3.4 + 3.0*fabs(x) + 2.0*fabs(y);
    } // time_delay

    static
    const char* time_units(void) {
        return "s";
    } // time_units

    static
    double point_force_3d_x(const double x,
                            const double y,
                            const double z) {
        return -0.3*x*z + 0.1*x*z;
    } // point_force_3d_x

    static
    double point_force_3d_y(const double x,
                            const double y,
                            const double z) {
        return -0.8*x*z + 0.2*y*y;
    } // point_force_3d_y

    static
    double point_force_3d_z(const double x,
                            const double y,
                            const double z) {
        return -0.3*x*x + 0.5*y*z;
    } // point_force_3d_z

    static
    double time_delay_3d(const double x,
                         const double y,
                         const double z) {
        return 3.4 + 3.0*fabs(x) + 2.0*fabs(z);
    } // time_delay

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestAuxiliaryFactoryPointForce::Tri::testAdd", "[TestAuxiliaryFactoryPointForce][add]") {
    pylith::sources::TestAuxiliaryFactoryPointForce(pylith::sources::TestAuxiliaryFactoryPointForce_Cases::Tri()).testAdd();
}
TEST_CASE("TestAuxiliaryFactoryPointForce::Tri::testSetValuesFromDB", "[TestAuxiliaryFactoryPointForce][testSetValuesFromDB]") {
    pylith::sources::TestAuxiliaryFactoryPointForce(pylith::sources::TestAuxiliaryFactoryPointForce_Cases::Tri()).testSetValuesFromDB();
}

TEST_CASE("TestAuxiliaryFactoryPointForce::Hex::testAdd", "[TestAuxiliaryFactoryPointForce][add]") {
    pylith::sources::TestAuxiliaryFactoryPointForce(pylith::sources::TestAuxiliaryFactoryPointForce_Cases::Hex()).testAdd();
}
TEST_CASE("TestAuxiliaryFactoryPointForce::Hex::testSetValuesFromDB", "[TestAuxiliaryFactoryPointForce][testSetValuesFromDB]") {
    pylith::sources::TestAuxiliaryFactoryPointForce(pylith::sources::TestAuxiliaryFactoryPointForce_Cases::Hex()).testSetValuesFromDB();
}

// --------------------------------------------------------------------------------------------------------------------
pylith::sources::TestAuxiliaryFactoryPointForce_Data*
pylith::sources::TestAuxiliaryFactoryPointForce_Cases::Tri(void) {
    pylith::sources::TestAuxiliaryFactoryPointForce_Data* data = new pylith::sources::TestAuxiliaryFactoryPointForce_Data();
    assert(data);

    data->auxDim = 2;
    data->dimension = 2;
    data->meshFilename = "data/tri_onecell.mesh";
    data->cs = new spatialdata::geocoords::CSCart();assert(data->cs);
    data->cs->setSpaceDim(data->dimension);

    assert(data->auxiliaryDB);
    data->auxiliaryDB->addValue("point_force_x", point_force_2d_x, point_force_units());
    data->auxiliaryDB->addValue("point_force_y", point_force_2d_y, point_force_units());
    data->auxiliaryDB->addValue("time_delay", time_delay_2d, time_units());
    data->auxiliaryDB->setDescription("auxiliary");
    data->auxiliaryDB->setCoordSys(*data->cs);

    return data;
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::sources::TestAuxiliaryFactoryPointForce_Data*
pylith::sources::TestAuxiliaryFactoryPointForce_Cases::Hex(void) {
    pylith::sources::TestAuxiliaryFactoryPointForce_Data* data = new pylith::sources::TestAuxiliaryFactoryPointForce_Data();
    assert(data);

    data->auxDim = 3;
    data->dimension = 3;
    data->meshFilename = "data/hex_onecell.mesh";
    data->cs = new spatialdata::geocoords::CSCart();assert(data->cs);
    data->cs->setSpaceDim(data->dimension);

    assert(data->auxiliaryDB);
    data->auxiliaryDB->addValue("point_force_x", point_force_3d_x, point_force_units());
    data->auxiliaryDB->addValue("point_force_y", point_force_3d_y, point_force_units());
    data->auxiliaryDB->addValue("point_force_z", point_force_3d_z, point_force_units());
    data->auxiliaryDB->addValue("time_delay", time_delay_3d, time_units());
    data->auxiliaryDB->setDescription("auxiliary");
    data->auxiliaryDB->setCoordSys(*data->cs);

    return data;
} // Hex
