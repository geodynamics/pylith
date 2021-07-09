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

#include "TestAuxiliaryFactoryElasticity.hh" // Implementation of cases

#include "pylith/materials/AuxiliaryFactoryElasticity.hh" // USES AuxiliaryFactoryElasticity
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cmath> // USES fabs()

// forward declarations
namespace pylith {
    namespace materials {
        class TestAuxiliaryFactoryElasticity_Tri;
        class TestAuxiliaryFactoryElasticity_Hex;
    } // materials
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::materials::TestAuxiliaryFactoryElasticity_Tri : public pylith::materials::TestAuxiliaryFactoryElasticity {
    CPPUNIT_TEST_SUB_SUITE(TestAuxiliaryFactoryElasticity_Tri, TestAuxiliaryFactoryElasticity);
    CPPUNIT_TEST_SUITE_END();

    static
    double density(const double x,
                   const double y) {
        return 3.4 + 3.0*fabs(x) + 2.0*fabs(y);
    } // density

    static
    const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    static
    double body_force_x(const double x,
                        const double y) {
        return -0.3*x*y + 0.2*y*y;
    } // body_force_x

    static
    double body_force_y(const double x,
                        const double y) {
        return +0.3*x*x + 0.2*x*y;
    } // body_force_y

    static
    const char* body_force_units(void) {
        return "kg/(m**2*s**2)";
    } // body_force_units

    static
    double gravity_field_x(const double x,
                           const double y) {
        return 0.0;
    } // gravity_field_x

    static
    double gravity_field_y(const double x,
                           const double y) {
        return -9.80665;
    } // gravity_field_y

    static
    const char* gravity_field_units(void) {
        return "m/s**2";
    } // gravity_field_units

protected:

    void setUp(void) {
        _auxDim = 2;
        TestAuxiliaryFactoryElasticity::setUp();

        CPPUNIT_ASSERT(_data);
        _data->dimension = 2;
        _data->meshFilename = "data/tri.mesh";
        _data->cs = new spatialdata::geocoords::CSCart();CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->dimension);

        _data->gravityField->setGravityDir(0.0, -1.0, 0.0);

        CPPUNIT_ASSERT(_data->auxiliaryDB);
        _data->auxiliaryDB->addValue("density", density, density_units());
        _data->auxiliaryDB->addValue("body_force_x", body_force_x, body_force_units());
        _data->auxiliaryDB->addValue("body_force_y", body_force_y, body_force_units());
        _data->auxiliaryDB->addValue("gravitational_acceleration_x", gravity_field_x, gravity_field_units());
        _data->auxiliaryDB->addValue("gravitational_acceleration_y", gravity_field_y, gravity_field_units());
        _data->auxiliaryDB->setLabel("auxiliary");
        _data->auxiliaryDB->setCoordSys(*_data->cs);

        _data->subfields["body_force"].description.numComponents = 2;
        _data->subfields["gravitational_acceleration"].description.numComponents = 2;

        _initialize();
    } // setUp

}; // TestAuxiliaryFactoryElasticity_Tri
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestAuxiliaryFactoryElasticity_Tri);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::materials::TestAuxiliaryFactoryElasticity_Hex : public pylith::materials::TestAuxiliaryFactoryElasticity {
    CPPUNIT_TEST_SUB_SUITE(TestAuxiliaryFactoryElasticity_Hex, TestAuxiliaryFactoryElasticity);
    CPPUNIT_TEST_SUITE_END();

    static
    double density(const double x,
                   const double y,
                   const double z) {
        return 6.4 + 3.0*fabs(x) + 2.0*fabs(y);
    } // density

    static
    const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    static
    double body_force_x(const double x,
                        const double y,
                        const double z) {
        return -0.3*x*y + 0.2*y*y;
    } // body_force_x

    static
    double body_force_y(const double x,
                        const double y,
                        const double z) {
        return +0.3*x*x + 0.2*x*y;
    } // body_force_y

    static
    double body_force_z(const double x,
                        const double y,
                        const double z) {
        return +0.3*x*y + 0.2*x*z;
    } // body_force_z

    static
    const char* body_force_units(void) {
        return "kg/(m**2*s**2)";
    } // body_force_units

    static
    double gravity_field_x(const double x,
                           const double y,
                           const double z) {
        return 0.0;
    } // gravity_field_x

    static
    double gravity_field_y(const double x,
                           const double y,
                           const double z) {
        return 0;
    } // gravity_field_y

    static
    double gravity_field_z(const double x,
                           const double y,
                           const double z) {
        return -9.80665;
    } // gravity_field_z

    static
    const char* gravity_field_units(void) {
        return "m/s**2";
    } // gravity_field_units

protected:

    void setUp(void) {
        _auxDim = 3;
        TestAuxiliaryFactoryElasticity::setUp();

        CPPUNIT_ASSERT(_data);
        _data->dimension = 3;
        _data->meshFilename = "data/hex.mesh";
        _data->cs = new spatialdata::geocoords::CSCart();CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->dimension);

        CPPUNIT_ASSERT(_data->auxiliaryDB);
        _data->auxiliaryDB->addValue("density", density, density_units());
        _data->auxiliaryDB->addValue("body_force_x", body_force_x, body_force_units());
        _data->auxiliaryDB->addValue("body_force_y", body_force_y, body_force_units());
        _data->auxiliaryDB->addValue("body_force_z", body_force_z, body_force_units());
        _data->auxiliaryDB->addValue("gravitational_acceleration_x", gravity_field_x, gravity_field_units());
        _data->auxiliaryDB->addValue("gravitational_acceleration_y", gravity_field_y, gravity_field_units());
        _data->auxiliaryDB->addValue("gravitational_acceleration_z", gravity_field_z, gravity_field_units());
        _data->auxiliaryDB->setLabel("auxiliary");
        _data->auxiliaryDB->setCoordSys(*_data->cs);

        _data->gravityField->setGravityDir(0.0, 0.0, -1.0);

        _data->subfields["body_force"].description.numComponents = 3;
        _data->subfields["gravitational_acceleration"].description.numComponents = 3;

        _initialize();
    } // setUp

}; // TestAuxiliaryFactoryElasticity_Hex
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestAuxiliaryFactoryElasticity_Hex);

// End of file
