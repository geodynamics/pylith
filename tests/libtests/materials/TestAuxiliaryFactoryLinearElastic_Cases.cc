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

#include "TestAuxiliaryFactoryLinearElastic.hh" // Implementation of cases

#include "pylith/materials/AuxiliaryFactoryElastic.hh" // USES AuxiliaryFactoryElastic
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "catch2/catch_test_macros.hpp"

#include <cmath> // USES fabs()

// forward declarations
namespace pylith {
    namespace materials {
        class TestAuxiliaryFactoryLinearElastic_Cases;
    } // materials
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::materials::TestAuxiliaryFactoryLinearElastic_Cases {
public:

    // Data factory methods
    static TestAuxiliaryFactoryLinearElastic_Data* Tri(void);

    static TestAuxiliaryFactoryLinearElastic_Data* Hex(void);

private:

    static
    double density_2d(const double x,
                      const double y) {
        return 3.4 + 3.0*fabs(x) + 2.0*fabs(y);
    } // density

    static
    const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    static
    double vs_2d(const double x,
                 const double y) {
        return 3.4 + 3.0*fabs(x) + 2.0*fabs(y);
    } // vs

    static
    const char* vs_units(void) {
        return "km/s";
    } // vs_units

    static
    double vp_2d(const double x,
                 const double y) {
        return 6.4 + 5.0*fabs(x) + 4.3*fabs(y);
    } // vp

    static
    const char* vp_units(void) {
        return "km/s";
    } // vp_units

    static
    double shear_modulus_2d(const double x,
                            const double y) {
        const double vsV = vs_2d(x,y);
        return density_2d(x,y) * (vsV*vsV);
    } // shear_modulus

    static
    double bulk_modulus_2d(const double x,
                           const double y) {
        const double vsV = vs_2d(x,y);
        const double vpV = vp_2d(x,y);
        return density_2d(x,y)*(vpV*vpV - 4.0/3.0*vsV*vsV);
    } // bulk_modulus

    static
    const char* modulus_units(void) {
        return "Pa";
    } // modulus_units

    static
    double reference_stress_2d_xx(const double x,
                                  const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // reference_stress_xx

    static
    double reference_stress_2d_yy(const double x,
                                  const double y) {
        return -0.8*x*y + 0.2*y*y;
    } // reference_stress_2d_yy

    static
    double reference_stress_2d_zz(const double x,
                                  const double y) {
        return -0.3*x*x + 0.5*y*y;
    } // reference_stress_2d_zz

    static
    double reference_stress_2d_xy(const double x,
                                  const double y) {
        return -0.3*x*x + 0.2*x*y;
    } // reference_stress_xy

    static
    const char* stress_units(void) {
        return "Pa";
    } // stress_units

    static
    double reference_strain_2d_xx(const double x,
                                  const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // reference_strain_2d_xx

    static
    double reference_strain_2d_yy(const double x,
                                  const double y) {
        return -0.8*x*y + 0.2*y*y;
    } // reference_strain_2d_yy

    static
    double reference_strain_2d_zz(const double x,
                                  const double y) {
        return -0.3*x*x + 0.5*y*y;
    } // reference_strain_2d_zz

    static
    double reference_strain_2d_xy(const double x,
                                  const double y) {
        return -0.3*x*x + 0.2*x*y;
    } // reference_strain_xy

    static
    const char* strain_units(void) {
        return "None";
    } // stress_units

    static
    double density_3d(const double x,
                      const double y,
                      const double z) {
        return 3.4 + 3.0*fabs(x) + 2.0*fabs(z);
    } // density

    static
    double vs_3d(const double x,
                 const double y,
                 const double z) {
        return 3.4 + 3.0*fabs(x) + 2.0*fabs(z);
    } // vs

    static
    double vp_3d(const double x,
                 const double y,
                 const double z) {
        return 6.4 + 4.0*fabs(x) + 5.3*fabs(z);
    } // vp

    static
    double shear_modulus_3d(const double x,
                            const double y,
                            const double z) {
        const double vsV = vs_3d(x,y,z);
        return density_3d(x,y,z) *(vsV*vsV);
    } // shear_modulus

    static
    double bulk_modulus_3d(const double x,
                           const double y,
                           const double z) {
        const double vsV = vs_3d(x,y,z);
        const double vpV = vp_3d(x,y,z);
        return density_3d(x,y,z)*(vpV*vpV - 4.0/3.0*vsV*vsV);
    } // bulk_modulus

    static
    double reference_stress_3d_xx(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*x*z + 0.1*x*z;
    } // reference_stress_3d_xx

    static
    double reference_stress_3d_yy(const double x,
                                  const double y,
                                  const double z) {
        return -0.8*x*z + 0.2*y*y;
    } // reference_stress_3d_yy

    static
    double reference_stress_3d_zz(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*x*x + 0.5*y*z;
    } // reference_stress_3d_zz

    static
    double reference_stress_3d_xy(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*x + 0.2*x*y;
    } // reference_stress_3d_xy

    static
    double reference_stress_3d_yz(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*x*z + 0.2*x*y;
    } // reference_stress_3d_yz

    static
    double reference_stress_3d_xz(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*x + 0.2*y;
    } // reference_stress_3d_xz

    static
    double reference_strain_3d_xx(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*x*x + 0.1*x*z;
    } // reference_strain_3d_xx

    static
    double reference_strain_3d_yy(const double x,
                                  const double y,
                                  const double z) {
        return -0.8*x*y + 0.2*y*y;
    } // reference_strain_3d_yy

    static
    double reference_strain_3d_zz(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*x*z + 0.5*y*y;
    } // reference_strain_3d_zz

    static
    double reference_strain_3d_xy(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*x*x + 0.2*x*z;
    } // reference_strain_3d_xy

    static
    double reference_strain_3d_yz(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*y*z + 0.2*x*z;
    } // reference_strain_3d_yz

    static
    double reference_strain_3d_xz(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*y*y + 0.2*x*z;
    } // reference_strain_3d_xz

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestAuxiliaryFactoryLinearElastic::Tri::testAdd", "[TestAuxiliaryFactoryLinearElastic][add]") {
    pylith::materials::TestAuxiliaryFactoryLinearElastic(pylith::materials::TestAuxiliaryFactoryLinearElastic_Cases::Tri()).testAdd();
}
TEST_CASE("TestAuxiliaryFactoryLinearElastic::Tri::testSetValuesFromDB", "[TestAuxiliaryFactoryLinearElastic][testSetValuesFromDB]") {
    pylith::materials::TestAuxiliaryFactoryLinearElastic(pylith::materials::TestAuxiliaryFactoryLinearElastic_Cases::Tri()).testSetValuesFromDB();
}

TEST_CASE("TestAuxiliaryFactoryLinearElastic::Hex::testAdd", "[TestAuxiliaryFactoryLinearElastic][add]") {
    pylith::materials::TestAuxiliaryFactoryLinearElastic(pylith::materials::TestAuxiliaryFactoryLinearElastic_Cases::Hex()).testAdd();
}
TEST_CASE("TestAuxiliaryFactoryLinearElastic::Hex::testSetValuesFromDB", "[TestAuxiliaryFactoryLinearElastic][testSetValuesFromDB]") {
    pylith::materials::TestAuxiliaryFactoryLinearElastic(pylith::materials::TestAuxiliaryFactoryLinearElastic_Cases::Hex()).testSetValuesFromDB();
}

// --------------------------------------------------------------------------------------------------------------------
pylith::materials::TestAuxiliaryFactoryLinearElastic_Data*
pylith::materials::TestAuxiliaryFactoryLinearElastic_Cases::Tri(void) {
    pylith::materials::TestAuxiliaryFactoryLinearElastic_Data* data = new pylith::materials::TestAuxiliaryFactoryLinearElastic_Data();
    assert(data);

    data->auxDim = 2;
    data->dimension = 2;
    data->meshFilename = "data/tri.mesh";
    data->cs = new spatialdata::geocoords::CSCart();assert(data->cs);
    data->cs->setSpaceDim(data->dimension);

    assert(data->auxiliaryDB);
    data->auxiliaryDB->addValue("density", density_2d, density_units());
    data->auxiliaryDB->addValue("vs", vs_2d, vs_units());
    data->auxiliaryDB->addValue("vp", vp_2d, vp_units());
    data->auxiliaryDB->addValue("shear_modulus", shear_modulus_2d, modulus_units());
    data->auxiliaryDB->addValue("bulk_modulus", bulk_modulus_2d, modulus_units());
    data->auxiliaryDB->addValue("reference_stress_xx", reference_stress_2d_xx, stress_units());
    data->auxiliaryDB->addValue("reference_stress_yy", reference_stress_2d_yy, stress_units());
    data->auxiliaryDB->addValue("reference_stress_zz", reference_stress_2d_zz, stress_units());
    data->auxiliaryDB->addValue("reference_stress_xy", reference_stress_2d_xy, stress_units());
    data->auxiliaryDB->addValue("reference_strain_xx", reference_strain_2d_xx, strain_units());
    data->auxiliaryDB->addValue("reference_strain_yy", reference_strain_2d_yy, strain_units());
    data->auxiliaryDB->addValue("reference_strain_zz", reference_strain_2d_zz, strain_units());
    data->auxiliaryDB->addValue("reference_strain_xy", reference_strain_2d_xy, strain_units());
    data->auxiliaryDB->setDescription("auxiliary");
    data->auxiliaryDB->setCoordSys(*data->cs);

    return data;
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::materials::TestAuxiliaryFactoryLinearElastic_Data*
pylith::materials::TestAuxiliaryFactoryLinearElastic_Cases::Hex(void) {
    pylith::materials::TestAuxiliaryFactoryLinearElastic_Data* data = new pylith::materials::TestAuxiliaryFactoryLinearElastic_Data();
    assert(data);

    data->auxDim = 3;
    data->dimension = 3;
    data->meshFilename = "data/hex.mesh";
    data->cs = new spatialdata::geocoords::CSCart();assert(data->cs);
    data->cs->setSpaceDim(data->dimension);

    assert(data->auxiliaryDB);
    data->auxiliaryDB->addValue("density", density_3d, density_units());
    data->auxiliaryDB->addValue("vs", vs_3d, vs_units());
    data->auxiliaryDB->addValue("vp", vp_3d, vp_units());
    data->auxiliaryDB->addValue("shear_modulus", shear_modulus_3d, modulus_units());
    data->auxiliaryDB->addValue("bulk_modulus", bulk_modulus_3d, modulus_units());
    data->auxiliaryDB->addValue("reference_stress_xx", reference_stress_3d_xx, stress_units());
    data->auxiliaryDB->addValue("reference_stress_yy", reference_stress_3d_yy, stress_units());
    data->auxiliaryDB->addValue("reference_stress_zz", reference_stress_3d_zz, stress_units());
    data->auxiliaryDB->addValue("reference_stress_xy", reference_stress_3d_xy, stress_units());
    data->auxiliaryDB->addValue("reference_stress_yz", reference_stress_3d_yz, stress_units());
    data->auxiliaryDB->addValue("reference_stress_xz", reference_stress_3d_xz, stress_units());
    data->auxiliaryDB->addValue("reference_strain_xx", reference_strain_3d_xx, strain_units());
    data->auxiliaryDB->addValue("reference_strain_yy", reference_strain_3d_yy, strain_units());
    data->auxiliaryDB->addValue("reference_strain_zz", reference_strain_3d_zz, strain_units());
    data->auxiliaryDB->addValue("reference_strain_xy", reference_strain_3d_xy, strain_units());
    data->auxiliaryDB->addValue("reference_strain_yz", reference_strain_3d_yz, strain_units());
    data->auxiliaryDB->addValue("reference_strain_xz", reference_strain_3d_xz, strain_units());
    data->auxiliaryDB->setDescription("auxiliary");
    data->auxiliaryDB->setCoordSys(*data->cs);

    return data;
} // Hex


// End of file
