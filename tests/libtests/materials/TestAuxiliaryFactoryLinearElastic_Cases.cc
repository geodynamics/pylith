// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestAuxiliaryFactoryLinearElastic.hh" // Implementation of cases

#include "pylith/materials/AuxiliaryFactoryElastic.hh" // USES AuxiliaryFactoryElastic
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Scales.hh" // USES Scales

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

    static const PylithReal LENGTH_SCALE;
    static const PylithReal KM;
    static const PylithReal TIME_SCALE;
    static const PylithReal PRESSURE_SCALE;
    static const PylithReal DENSITY_SCALE;

private:

    static
    double density_2d(const double x,
                      const double y) {
        return 2500.0 + 3.0*fabs(x)/KM + 2.0*fabs(y)/KM;
    } // density

    static
    const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    static
    double vs_2d(const double x,
                 const double y) {
        return 1000.0 + 30.0*fabs(x)/KM + 20.0*fabs(y)/KM;
    } // vs

    static
    const char* vs_units(void) {
        return "m/s";
    } // vs_units

    static
    double vp_2d(const double x,
                 const double y) {
        return 2000.0 + 50.0*fabs(x)/KM + 43.0*fabs(y)/KM;
    } // vp

    static
    const char* vp_units(void) {
        return "m/s";
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
        return -0.3*x*x/(KM*KM) + 0.1*x*y/(KM*KM);
    } // reference_stress_xx

    static
    double reference_stress_2d_yy(const double x,
                                  const double y) {
        return -8.0e+6*x*y/(KM*KM) + 2.0e+6*y*y/(KM*KM);
    } // reference_stress_2d_yy

    static
    double reference_stress_2d_zz(const double x,
                                  const double y) {
        return -3.0e+6*x*x/(KM*KM) + 5.0e+6*y*y/(KM*KM);
    } // reference_stress_2d_zz

    static
    double reference_stress_2d_xy(const double x,
                                  const double y) {
        return -3.0e+6*x*x/(KM*KM) + 2.0e+6*x*y/(KM*KM);
    } // reference_stress_xy

    static
    const char* stress_units(void) {
        return "Pa";
    } // stress_units

    static
    double reference_strain_2d_xx(const double x,
                                  const double y) {
        return -0.3*x*x/(KM*KM) + 0.1*x*y/(KM*KM);
    } // reference_strain_2d_xx

    static
    double reference_strain_2d_yy(const double x,
                                  const double y) {
        return -0.8*x*y/(KM*KM) + 0.2*y*y/(KM*KM);
    } // reference_strain_2d_yy

    static
    double reference_strain_2d_zz(const double x,
                                  const double y) {
        return -0.3*x*x/(KM*KM) + 0.5*y*y/(KM*KM);
    } // reference_strain_2d_zz

    static
    double reference_strain_2d_xy(const double x,
                                  const double y) {
        return -0.3*x*x/(KM*KM) + 0.2*x*y/(KM*KM);
    } // reference_strain_xy

    static
    const char* strain_units(void) {
        return "None";
    } // stress_units

    static
    double density_3d(const double x,
                      const double y,
                      const double z) {
        return 2500.0 + 3.0*fabs(x)/KM + 2.0*fabs(z)/KM;
    } // density

    static
    double vs_3d(const double x,
                 const double y,
                 const double z) {
        return 1000.0 + 300.0*fabs(x)/KM + 200.0*fabs(z)/KM;
    } // vs

    static
    double vp_3d(const double x,
                 const double y,
                 const double z) {
        return 2500.0 + 400.0*fabs(x)/KM + 5.3*fabs(z)/KM;
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
        return -3.0e+6*x*z/(KM*KM) + 1.0e+6*x*z/(KM*KM);
    } // reference_stress_3d_xx

    static
    double reference_stress_3d_yy(const double x,
                                  const double y,
                                  const double z) {
        return -8.0e+6*x*z/(KM*KM) + 2.0e+6*y*y/(KM*KM);
    } // reference_stress_3d_yy

    static
    double reference_stress_3d_zz(const double x,
                                  const double y,
                                  const double z) {
        return -3.0e+6*x*x/(KM*KM) + 5.0e+6*y*z/(KM*KM);
    } // reference_stress_3d_zz

    static
    double reference_stress_3d_xy(const double x,
                                  const double y,
                                  const double z) {
        return -3.0e+6*x/KM + 2.0e+6*x*y/(KM*KM);
    } // reference_stress_3d_xy

    static
    double reference_stress_3d_yz(const double x,
                                  const double y,
                                  const double z) {
        return -3.0e+6*x*z/(KM*KM) + 2.0e+6*x*y/(KM*KM);
    } // reference_stress_3d_yz

    static
    double reference_stress_3d_xz(const double x,
                                  const double y,
                                  const double z) {
        return -3.0e+6*x/KM + 2.0e+6*y/KM;
    } // reference_stress_3d_xz

    static
    double reference_strain_3d_xx(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*x*x/(KM*KM) + 0.1*x*z/(KM*KM);
    } // reference_strain_3d_xx

    static
    double reference_strain_3d_yy(const double x,
                                  const double y,
                                  const double z) {
        return -0.8*x*y/(KM*KM) + 0.2*y*y/(KM*KM);
    } // reference_strain_3d_yy

    static
    double reference_strain_3d_zz(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*x*z/(KM*KM) + 0.5*y*y/(KM*KM);
    } // reference_strain_3d_zz

    static
    double reference_strain_3d_xy(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*x*x/(KM*KM) + 0.2*x*z/(KM*KM);
    } // reference_strain_3d_xy

    static
    double reference_strain_3d_yz(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*y*z/(KM*KM) + 0.2*x*z/(KM*KM);
    } // reference_strain_3d_yz

    static
    double reference_strain_3d_xz(const double x,
                                  const double y,
                                  const double z) {
        return -0.3*y*y/(KM*KM) + 0.2*x*z/(KM*KM);
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

const PylithReal pylith::materials::TestAuxiliaryFactoryLinearElastic_Cases::LENGTH_SCALE = 1.0;
const PylithReal pylith::materials::TestAuxiliaryFactoryLinearElastic_Cases::TIME_SCALE = 2.0;
const PylithReal pylith::materials::TestAuxiliaryFactoryLinearElastic_Cases::PRESSURE_SCALE = 3.0e+6;
const PylithReal pylith::materials::TestAuxiliaryFactoryLinearElastic_Cases::KM = 1.0e+3;

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

    assert(data->scales);
    data->scales->setLengthScale(LENGTH_SCALE);
    data->scales->setTimeScale(TIME_SCALE);
    data->scales->setPressureScale(PRESSURE_SCALE);

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

    assert(data->scales);
    data->scales->setLengthScale(LENGTH_SCALE);
    data->scales->setTimeScale(TIME_SCALE);
    data->scales->setPressureScale(PRESSURE_SCALE);

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
