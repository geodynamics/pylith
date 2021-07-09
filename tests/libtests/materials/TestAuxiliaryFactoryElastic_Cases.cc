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

#include "TestAuxiliaryFactoryElastic.hh" // Implementation of cases

#include "pylith/materials/AuxiliaryFactoryElastic.hh" // USES AuxiliaryFactoryElastic
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cmath> // USES fabs()

// forward declarations
namespace pylith {
    namespace materials {
        class TestAuxiliaryFactoryElastic_Tri;
        class TestAuxiliaryFactoryElastic_Hex;
    } // materials
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::materials::TestAuxiliaryFactoryElastic_Tri : public pylith::materials::TestAuxiliaryFactoryElastic {
    CPPUNIT_TEST_SUB_SUITE(TestAuxiliaryFactoryElastic_Tri, TestAuxiliaryFactoryElastic);
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
    double vs(const double x,
              const double y) {
        return 3.4 + 3.0*fabs(x) + 2.0*fabs(y);
    } // vs

    static
    const char* vs_units(void) {
        return "km/s";
    } // vs_units

    static
    double vp(const double x,
              const double y) {
        return 6.4 + 5.0*fabs(x) + 4.3*fabs(y);
    } // vp

    static
    const char* vp_units(void) {
        return "km/s";
    } // vp_units

    static
    double shear_modulus(const double x,
                         const double y) {
        const double vsV = vs(x,y);
        return density(x,y) * (vsV*vsV);
    } // shear_modulus

    static
    double bulk_modulus(const double x,
                        const double y) {
        const double vsV = vs(x,y);
        const double vpV = vp(x,y);
        return density(x,y)*(vpV*vpV - 4.0/3.0*vsV*vsV);
    } // bulk_modulus

    static
    const char* modulus_units(void) {
        return "Pa";
    } // modulus_units

    static
    double reference_stress_xx(const double x,
                               const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // reference_stress_xx

    static
    double reference_stress_yy(const double x,
                               const double y) {
        return -0.8*x*y + 0.2*y*y;
    } // reference_stress_yy

    static
    double reference_stress_zz(const double x,
                               const double y) {
        return -0.3*x*x + 0.5*y*y;
    } // reference_stress_zz

    static
    double reference_stress_xy(const double x,
                               const double y) {
        return -0.3*x*x + 0.2*x*y;
    } // reference_stress_xy

    static
    const char* stress_units(void) {
        return "Pa";
    } // stress_units

    static
    double reference_strain_xx(const double x,
                               const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // reference_strain_xx

    static
    double reference_strain_yy(const double x,
                               const double y) {
        return -0.8*x*y + 0.2*y*y;
    } // reference_strain_yy

    static
    double reference_strain_zz(const double x,
                               const double y) {
        return -0.3*x*x + 0.5*y*y;
    } // reference_strain_zz

    static
    double reference_strain_xy(const double x,
                               const double y) {
        return -0.3*x*x + 0.2*x*y;
    } // reference_strain_xy

    static
    const char* strain_units(void) {
        return "None";
    } // stress_units

protected:

    void setUp(void) {
        _auxDim = 2;
        TestAuxiliaryFactoryElastic::setUp();

        CPPUNIT_ASSERT(_data);
        _data->dimension = 2;
        _data->meshFilename = "data/tri.mesh";
        _data->cs = new spatialdata::geocoords::CSCart();CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->dimension);

        CPPUNIT_ASSERT(_data->auxiliaryDB);
        _data->auxiliaryDB->addValue("density", density, density_units());
        _data->auxiliaryDB->addValue("vs", vs, vs_units());
        _data->auxiliaryDB->addValue("vp", vp, vp_units());
        _data->auxiliaryDB->addValue("shear_modulus", shear_modulus, modulus_units());
        _data->auxiliaryDB->addValue("bulk_modulus", bulk_modulus, modulus_units());
        _data->auxiliaryDB->addValue("reference_stress_xx", reference_stress_xx, stress_units());
        _data->auxiliaryDB->addValue("reference_stress_yy", reference_stress_xx, stress_units());
        _data->auxiliaryDB->addValue("reference_stress_zz", reference_stress_xx, stress_units());
        _data->auxiliaryDB->addValue("reference_stress_xy", reference_stress_xx, stress_units());
        _data->auxiliaryDB->addValue("reference_strain_xx", reference_strain_xx, strain_units());
        _data->auxiliaryDB->addValue("reference_strain_yy", reference_strain_xx, strain_units());
        _data->auxiliaryDB->addValue("reference_strain_zz", reference_strain_xx, strain_units());
        _data->auxiliaryDB->addValue("reference_strain_xy", reference_strain_xx, strain_units());
        _data->auxiliaryDB->setLabel("auxiliary");
        _data->auxiliaryDB->setCoordSys(*_data->cs);

        _data->subfields["reference_stress"].description.numComponents = 4;
        _data->subfields["reference_stress"].description.vectorFieldType = pylith::topology::Field::OTHER;
        _data->subfields["reference_strain"].description.numComponents = 4;
        _data->subfields["reference_strain"].description.vectorFieldType = pylith::topology::Field::OTHER;

        _initialize();
    } // setUp

}; // TestAuxiliaryFactoryElastic_Tri
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestAuxiliaryFactoryElastic_Tri);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::materials::TestAuxiliaryFactoryElastic_Hex : public pylith::materials::TestAuxiliaryFactoryElastic {
    CPPUNIT_TEST_SUB_SUITE(TestAuxiliaryFactoryElastic_Hex, TestAuxiliaryFactoryElastic);
    CPPUNIT_TEST_SUITE_END();

    static
    double density(const double x,
                   const double y,
                   const double z) {
        return 3.4 + 3.0*fabs(x) + 2.0*fabs(z);
    } // density

    static
    const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    static
    double vs(const double x,
              const double y,
              const double z) {
        return 3.4 + 3.0*fabs(x) + 2.0*fabs(z);
    } // vs

    static
    const char* vs_units(void) {
        return "km/s";
    } // vs_units

    static
    double vp(const double x,
              const double y,
              const double z) {
        return 6.4 + 4.0*fabs(x) + 5.3*fabs(z);
    } // vp

    static
    const char* vp_units(void) {
        return "km/s";
    } // vp_units

    static
    double shear_modulus(const double x,
                         const double y,
                         const double z) {
        const double vsV = vs(x,y,z);
        return density(x,y,z) *(vsV*vsV);
    } // shear_modulus

    static
    double bulk_modulus(const double x,
                        const double y,
                        const double z) {
        const double vsV = vs(x,y,z);
        const double vpV = vp(x,y,z);
        return density(x,y,z)*(vpV*vpV - 4.0/3.0*vsV*vsV);
    } // bulk_modulus

    static
    const char* modulus_units(void) {
        return "Pa";
    } // modulus_units

    static
    double reference_stress_xx(const double x,
                               const double y,
                               const double z) {
        return -0.3*x*z + 0.1*x*z;
    } // reference_stress_xx

    static
    double reference_stress_yy(const double x,
                               const double y,
                               const double z) {
        return -0.8*x*z + 0.2*y*y;
    } // reference_stress_yy

    static
    double reference_stress_zz(const double x,
                               const double y,
                               const double z) {
        return -0.3*x*x + 0.5*y*z;
    } // reference_stress_zz

    static
    double reference_stress_xy(const double x,
                               const double y,
                               const double z) {
        return -0.3*x + 0.2*x*y;
    } // reference_stress_xy

    static
    double reference_stress_yz(const double x,
                               const double y,
                               const double z) {
        return -0.3*x*z + 0.2*x*y;
    } // reference_stress_yz

    static
    double reference_stress_xz(const double x,
                               const double y,
                               const double z) {
        return -0.3*x + 0.2*y;
    } // reference_stress_xz

    static
    const char* stress_units(void) {
        return "Pa";
    } // stress_units

    static
    double reference_strain_xx(const double x,
                               const double y,
                               const double z) {
        return -0.3*x*x + 0.1*x*z;
    } // reference_strain_xx

    static
    double reference_strain_yy(const double x,
                               const double y,
                               const double z) {
        return -0.8*x*y + 0.2*y*y;
    } // reference_strain_yy

    static
    double reference_strain_zz(const double x,
                               const double y,
                               const double z) {
        return -0.3*x*z + 0.5*y*y;
    } // reference_strain_zz

    static
    double reference_strain_xy(const double x,
                               const double y,
                               const double z) {
        return -0.3*x*x + 0.2*x*z;
    } // reference_strain_xy

    static
    double reference_strain_yz(const double x,
                               const double y,
                               const double z) {
        return -0.3*y*z + 0.2*x*z;
    } // reference_strain_yz

    static
    double reference_strain_xz(const double x,
                               const double y,
                               const double z) {
        return -0.3*y*y + 0.2*x*z;
    } // reference_strain_xz

    static
    const char* strain_units(void) {
        return "None";
    } // stress_units

protected:

    void setUp(void) {
        _auxDim = 3;
        TestAuxiliaryFactoryElastic::setUp();

        CPPUNIT_ASSERT(_data);
        _data->dimension = 3;
        _data->meshFilename = "data/hex.mesh";
        _data->cs = new spatialdata::geocoords::CSCart();CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->dimension);

        CPPUNIT_ASSERT(_data->auxiliaryDB);
        _data->auxiliaryDB->addValue("density", density, density_units());
        _data->auxiliaryDB->addValue("vs", vs, vs_units());
        _data->auxiliaryDB->addValue("vp", vp, vp_units());
        _data->auxiliaryDB->addValue("shear_modulus", shear_modulus, modulus_units());
        _data->auxiliaryDB->addValue("bulk_modulus", bulk_modulus, modulus_units());
        _data->auxiliaryDB->addValue("reference_stress_xx", reference_stress_xx, stress_units());
        _data->auxiliaryDB->addValue("reference_stress_yy", reference_stress_xx, stress_units());
        _data->auxiliaryDB->addValue("reference_stress_zz", reference_stress_xx, stress_units());
        _data->auxiliaryDB->addValue("reference_stress_xy", reference_stress_xx, stress_units());
        _data->auxiliaryDB->addValue("reference_stress_yz", reference_stress_xx, stress_units());
        _data->auxiliaryDB->addValue("reference_stress_xz", reference_stress_xx, stress_units());
        _data->auxiliaryDB->addValue("reference_strain_xx", reference_strain_xx, strain_units());
        _data->auxiliaryDB->addValue("reference_strain_yy", reference_strain_xx, strain_units());
        _data->auxiliaryDB->addValue("reference_strain_zz", reference_strain_xx, strain_units());
        _data->auxiliaryDB->addValue("reference_strain_xy", reference_strain_xx, strain_units());
        _data->auxiliaryDB->addValue("reference_strain_yz", reference_strain_xx, strain_units());
        _data->auxiliaryDB->addValue("reference_strain_xz", reference_strain_xx, strain_units());
        _data->auxiliaryDB->setLabel("auxiliary");
        _data->auxiliaryDB->setCoordSys(*_data->cs);

        _initialize();
    } // setUp

}; // TestAuxiliaryFactoryElastic_Hex
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestAuxiliaryFactoryElastic_Hex);

// End of file
