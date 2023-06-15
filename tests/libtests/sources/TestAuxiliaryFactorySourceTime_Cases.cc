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

#include "TestAuxiliaryFactorySourceTime.hh" // Implementation of cases

#include "pylith/sources/AuxiliaryFactorySourceTime.hh" // USES AuxiliaryFactorySourceTime
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "catch2/catch_test_macros.hpp"

#include <cmath> // USES fabs()

// forward declarations
namespace pylith {
    namespace sources {
        class TestAuxiliaryFactorySourceTime_Cases;
    } // sources
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::sources::TestAuxiliaryFactorySourceTime_Cases {
public:

    // Data factory methods
    static TestAuxiliaryFactorySourceTime_Data* Tri(void);

    static TestAuxiliaryFactorySourceTime_Data* Hex(void);

private:

    static
    double moment_tensor_2d_xx(const double x,
                               const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // moment_tensor_2d_xx

    static
    double moment_tensor_2d_yy(const double x,
                               const double y) {
        return -0.8*x*y + 0.2*y*y;
    } // moment_tensor_2d_yy

    static
    double moment_tensor_2d_zz(const double x,
                               const double y) {
        return -0.3*x*x + 0.5*y*y;
    } // moment_tensor_2d_zz

    static
    double moment_tensor_2d_xy(const double x,
                               const double y) {
        return -0.3*x*x + 0.3*x*y;
    }

    static
    const char* pressure_units(void) {
        return "Pa";
    } // pressure_units

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
    double center_frequency_2d(const double x,
                               const double y) {
        return 6.4 + 5.0*fabs(x) + 4.3*fabs(y);
    } // center_frequency

    static
    const char* frequency_units(void) {
        return "1/s";
    } // frequency_units

    static
    double moment_tensor_3d_xx(const double x,
                               const double y,
                               const double z) {
        return -0.3*x*z + 0.1*x*z;
    } // moment_tensor_3d_xx

    static
    double moment_tensor_3d_yy(const double x,
                               const double y,
                               const double z) {
        return -0.8*x*z + 0.2*y*y;
    } // moment_tensor_3d_yy

    static
    double moment_tensor_3d_zz(const double x,
                               const double y,
                               const double z) {
        return -0.3*x*x + 0.5*y*z;
    } // moment_tensor_3d_zz

    static
    double moment_tensor_3d_xy(const double x,
                               const double y,
                               const double z) {
        return -0.3*x + 0.2*x*y;
    } // moment_tensor_3d_xy

    static
    double moment_tensor_3d_yz(const double x,
                               const double y,
                               const double z) {
        return -0.3*x*z + 0.2*x*y;
    } // moment_tensor_3d_yz

    static
    double moment_tensor_3d_xz(const double x,
                               const double y,
                               const double z) {
        return -0.3*x + 0.2*y;
    } // moment_tensor_3d_xz

    static
    double time_delay_3d(const double x,
                         const double y,
                         const double z) {
        return 3.4 + 3.0*fabs(x) + 2.0*fabs(z);
    } // time_delay

    static
    double center_frequency_3d(const double x,
                               const double y,
                               const double z) {
        return 6.4 + 4.0*fabs(x) + 5.3*fabs(z);
    } // center_frequency

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestAuxiliaryFactorySourceTime::Tri::testAdd", "[TestAuxiliaryFactorySourceTime][add]") {
    pylith::sources::TestAuxiliaryFactorySourceTime(pylith::sources::TestAuxiliaryFactorySourceTime_Cases::Tri()).testAdd();
}
TEST_CASE("TestAuxiliaryFactorySourceTime::Tri::testSetValuesFromDB", "[TestAuxiliaryFactorySourceTime][testSetValuesFromDB]") {
    pylith::sources::TestAuxiliaryFactorySourceTime(pylith::sources::TestAuxiliaryFactorySourceTime_Cases::Tri()).testSetValuesFromDB();
}

TEST_CASE("TestAuxiliaryFactorySourceTime::Hex::testAdd", "[TestAuxiliaryFactorySourceTime][add]") {
    pylith::sources::TestAuxiliaryFactorySourceTime(pylith::sources::TestAuxiliaryFactorySourceTime_Cases::Hex()).testAdd();
}
TEST_CASE("TestAuxiliaryFactorySourceTime::Hex::testSetValuesFromDB", "[TestAuxiliaryFactorySourceTime][testSetValuesFromDB]") {
    pylith::sources::TestAuxiliaryFactorySourceTime(pylith::sources::TestAuxiliaryFactorySourceTime_Cases::Hex()).testSetValuesFromDB();
}

// --------------------------------------------------------------------------------------------------------------------
pylith::sources::TestAuxiliaryFactorySourceTime_Data*
pylith::sources::TestAuxiliaryFactorySourceTime_Cases::Tri(void) {
    pylith::sources::TestAuxiliaryFactorySourceTime_Data* data = new pylith::sources::TestAuxiliaryFactorySourceTime_Data();
    assert(data);

    data->auxDim = 2;
    data->dimension = 2;
    data->meshFilename = "data/tri_onecell.mesh";
    data->cs = new spatialdata::geocoords::CSCart();assert(data->cs);
    data->cs->setSpaceDim(data->dimension);

    assert(data->auxiliaryDB);
    data->auxiliaryDB->addValue("moment_tensor_xx", moment_tensor_2d_xx, pressure_units());
    data->auxiliaryDB->addValue("moment_tensor_yy", moment_tensor_2d_yy, pressure_units());
    data->auxiliaryDB->addValue("moment_tensor_zz", moment_tensor_2d_zz, pressure_units());
    data->auxiliaryDB->addValue("moment_tensor_xy", moment_tensor_2d_xy, pressure_units());
    data->auxiliaryDB->addValue("time_delay", time_delay_2d, time_units());
    data->auxiliaryDB->addValue("center_frequency", center_frequency_2d, frequency_units());
    data->auxiliaryDB->setDescription("auxiliary");
    data->auxiliaryDB->setCoordSys(*data->cs);

    return data;
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::sources::TestAuxiliaryFactorySourceTime_Data*
pylith::sources::TestAuxiliaryFactorySourceTime_Cases::Hex(void) {
    pylith::sources::TestAuxiliaryFactorySourceTime_Data* data = new pylith::sources::TestAuxiliaryFactorySourceTime_Data();
    assert(data);

    data->auxDim = 3;
    data->dimension = 3;
    data->meshFilename = "data/hex_onecell.mesh";
    data->cs = new spatialdata::geocoords::CSCart();assert(data->cs);
    data->cs->setSpaceDim(data->dimension);

    assert(data->auxiliaryDB);
    data->auxiliaryDB->addValue("moment_tensor_xx", moment_tensor_3d_xx, pressure_units());
    data->auxiliaryDB->addValue("moment_tensor_yy", moment_tensor_3d_yy, pressure_units());
    data->auxiliaryDB->addValue("moment_tensor_zz", moment_tensor_3d_zz, pressure_units());
    data->auxiliaryDB->addValue("moment_tensor_xy", moment_tensor_3d_xy, pressure_units());
    data->auxiliaryDB->addValue("moment_tensor_yz", moment_tensor_3d_yz, pressure_units());
    data->auxiliaryDB->addValue("moment_tensor_xz", moment_tensor_3d_xz, pressure_units());
    data->auxiliaryDB->addValue("time_delay", time_delay_3d, time_units());
    data->auxiliaryDB->addValue("center_frequency", center_frequency_3d, frequency_units());
    data->auxiliaryDB->setDescription("auxiliary");
    data->auxiliaryDB->setCoordSys(*data->cs);

    return data;
} // Hex
