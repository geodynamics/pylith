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

#include "TestAuxiliaryFactoryElasticity.hh" // Implementation of cases

#include "pylith/materials/AuxiliaryFactoryElasticity.hh" // USES AuxiliaryFactoryElasticity
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "catch2/catch_test_macros.hpp"

#include <cmath> // USES fabs()

// forward declarations
namespace pylith {
    namespace materials {
        class TestAuxiliaryFactoryElasticity_Cases;
    } // materials
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::materials::TestAuxiliaryFactoryElasticity_Cases {
public:

    // Data factory methods
    static TestAuxiliaryFactoryElasticity_Data* Tri(void);

    static TestAuxiliaryFactoryElasticity_Data* Hex(void);

private:

    static
    double density_2d(const double x,
                      const double y) {
        return 6.4 + 3.0*fabs(x) + 2.0*fabs(y);
    } // density

    static
    double density_3d(const double x,
                      const double y,
                      const double z) {
        return 6.4 + 3.0*fabs(x) + 2.0*fabs(y) + 1.1*fabs(z);
    } // density

    static
    const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    static
    double body_force_2d_x(const double x,
                           const double y) {
        return -0.3*x*y + 0.2*y*y;
    } // body_force_x

    static
    double body_force_2d_y(const double x,
                           const double y) {
        return +0.3*x*x + 0.2*x*y;
    } // body_force_y

    static
    double body_force_3d_x(const double x,
                           const double y,
                           const double z) {
        return -0.3*x*y + 0.2*y*y;
    } // body_force_x

    static
    double body_force_3d_y(const double x,
                           const double y,
                           const double z) {
        return +0.3*x*x + 0.2*x*y;
    } // body_force_y

    static
    double body_force_3d_z(const double x,
                           const double y,
                           const double z) {
        return +0.3*x*y + 0.2*x*z;
    } // body_force_z

    static
    const char* body_force_units(void) {
        return "kg/(m**2*s**2)";
    } // body_force_units

    static
    double gravity_field_2d_x(const double x,
                              const double y) {
        return 0.0;
    } // gravity_field_x

    static
    double gravity_field_2d_y(const double x,
                              const double y) {
        return -9.80665;
    } // gravity_field_y

    static
    double gravity_field_3d_x(const double x,
                              const double y,
                              const double z) {
        return 0.0;
    } // gravity_field_x

    static
    double gravity_field_3d_y(const double x,
                              const double y,
                              const double z) {
        return 0.0;
    } // gravity_field_y

    static
    double gravity_field_3d_z(const double x,
                              const double y,
                              const double z) {
        return -9.80665;
    } // gravity_field_z

    static
    const char* gravity_field_units(void) {
        return "m/s**2";
    } // gravity_field_units

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestAuxiliaryFactoryElasticity::Tri::testAdd", "[TestAuxiliaryFactoryElasticity][add]") {
    pylith::materials::TestAuxiliaryFactoryElasticity(pylith::materials::TestAuxiliaryFactoryElasticity_Cases::Tri()).testAdd();
}
TEST_CASE("TestAuxiliaryFactoryElasticity::Tri::testSetValuesFromDB", "[TestAuxiliaryFactoryElasticity][testSetValuesFromDB]") {
    pylith::materials::TestAuxiliaryFactoryElasticity(pylith::materials::TestAuxiliaryFactoryElasticity_Cases::Tri()).testSetValuesFromDB();
}

TEST_CASE("TestAuxiliaryFactoryElasticity::Hex::testAdd", "[TestAuxiliaryFactoryElasticity][add]") {
    pylith::materials::TestAuxiliaryFactoryElasticity(pylith::materials::TestAuxiliaryFactoryElasticity_Cases::Hex()).testAdd();
}
TEST_CASE("TestAuxiliaryFactoryElasticity::Hex::testSetValuesFromDB", "[TestAuxiliaryFactoryElasticity][testSetValuesFromDB]") {
    pylith::materials::TestAuxiliaryFactoryElasticity(pylith::materials::TestAuxiliaryFactoryElasticity_Cases::Hex()).testSetValuesFromDB();
}

// ------------------------------------------------------------------------------------------------
pylith::materials::TestAuxiliaryFactoryElasticity_Data*
pylith::materials::TestAuxiliaryFactoryElasticity_Cases::Tri(void) {
    pylith::materials::TestAuxiliaryFactoryElasticity_Data* data = new pylith::materials::TestAuxiliaryFactoryElasticity_Data();
    assert(data);

    data->auxDim = 2;
    data->dimension = 2;
    data->meshFilename = "data/tri.mesh";
    data->cs = new spatialdata::geocoords::CSCart();assert(data->cs);
    data->cs->setSpaceDim(data->dimension);

    data->gravityField->setGravityDir(0.0, -1.0, 0.0);

    assert(data->auxiliaryDB);
    data->auxiliaryDB->addValue("density", density_2d, density_units());
    data->auxiliaryDB->addValue("body_force_x", body_force_2d_x, body_force_units());
    data->auxiliaryDB->addValue("body_force_y", body_force_2d_y, body_force_units());
    data->auxiliaryDB->addValue("gravitational_acceleration_x", gravity_field_2d_x, gravity_field_units());
    data->auxiliaryDB->addValue("gravitational_acceleration_y", gravity_field_2d_y, gravity_field_units());
    data->auxiliaryDB->setDescription("auxiliary");
    data->auxiliaryDB->setCoordSys(*data->cs);

    return data;
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::materials::TestAuxiliaryFactoryElasticity_Data*
pylith::materials::TestAuxiliaryFactoryElasticity_Cases::Hex(void) {
    pylith::materials::TestAuxiliaryFactoryElasticity_Data* data = new pylith::materials::TestAuxiliaryFactoryElasticity_Data();
    assert(data);

    data->auxDim = 3;
    data->dimension = 3;
    data->meshFilename = "data/hex.mesh";
    data->cs = new spatialdata::geocoords::CSCart();assert(data->cs);
    data->cs->setSpaceDim(data->dimension);

    assert(data->auxiliaryDB);
    data->auxiliaryDB->addValue("density", density_3d, density_units());
    data->auxiliaryDB->addValue("body_force_x", body_force_3d_x, body_force_units());
    data->auxiliaryDB->addValue("body_force_y", body_force_3d_y, body_force_units());
    data->auxiliaryDB->addValue("body_force_z", body_force_3d_z, body_force_units());
    data->auxiliaryDB->addValue("gravitational_acceleration_x", gravity_field_3d_x, gravity_field_units());
    data->auxiliaryDB->addValue("gravitational_acceleration_y", gravity_field_3d_y, gravity_field_units());
    data->auxiliaryDB->addValue("gravitational_acceleration_z", gravity_field_3d_z, gravity_field_units());
    data->auxiliaryDB->setDescription("auxiliary");
    data->auxiliaryDB->setCoordSys(*data->cs);

    data->gravityField->setGravityDir(0.0, 0.0, -1.0);

    data->subfields["body_force"].description.numComponents = 3;
    data->subfields["gravitational_acceleration"].description.numComponents = 3;

    return data;
} // Hex


// End of file
