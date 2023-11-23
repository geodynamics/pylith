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
    double fluid_density_2d(const double x,
                            const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // fluid_density_2d

    static
    const char* fluid_density_units(void) {
        return "kg/m**3";
    } // fluid_density_units

    static
    double fluid_viscosity_2d(const double x,
                              const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // fluid_viscosity_2d

    static
    const char* fluid_viscosity_units(void) {
        return "Pa*s";
    } // fluid_viscosity_units

    static
    double isotropic_permeability_2d(const double x,
                                     const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // isotropic_permeability_2d

    static
    const char* isotropic_permeability_units(void) {
        return "m**2";
    } // isotropic_permeability_units

    static
    double wellbore_radius_2d(const double x,
                              const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // wellbore_radius_2d

    static
    const char* wellbore_radius_units(void) {
        return "m";
    } // wellbore_radius_units

    static
    double wellbore_length_2d(const double x,
                              const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // wellbore_length_2d

    static
    const char* wellbore_length_units(void) {
        return "m";
    } // wellbore_length_units

    static
    double wellbore_pressure_2d(const double x,
                                const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // wellbore_pressure_2d

    static
    const char* wellbore_pressure_units(void) {
        return "Pa";
    } // wellbore_pressure_units

    static
    double wellbore_character_2d(const double x,
                                 const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // wellbore_character_2d

    static
    const char* wellbore_character_units(void) {
        return "one";
    } // wellbore_character_units

    static
    double time_delay_2d(const double x,
                         const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // time_delay_2d

    static
    const char* time_delay_units(void) {
        return "s";
    } // time_delay_units

    static
    double fluid_density_3d(const double x,
                            const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // fluid_density_3d

    static
    double fluid_viscosity_3d(const double x,
                              const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // fluid_viscosity_3d

    static
    double isotropic_permeability_3d(const double x,
                                     const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // isotropic_permeability_3d

    static
    double wellbore_radius_3d(const double x,
                              const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // wellbore_radius_3d

    static
    double wellbore_length_3d(const double x,
                              const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // wellbore_length_3d

    static
    double wellbore_pressure_3d(const double x,
                                const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // wellbore_pressure_3d

    static
    double wellbore_character_3d(const double x,
                                 const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // wellbore_character_2d

    static
    double time_delay_3d(const double x,
                         const double y) {
        return -0.3*x*x + 0.1*x*y;
    } // time_delay_3d

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
    data->auxiliaryDB->addValue("fluid_density", fluid_density_2d, fluid_density_units());
    data->auxiliaryDB->addValue("fluid_viscosity", fluid_viscosity_2d, fluid_viscosity_units());
    data->auxiliaryDB->addValue("isotropic_permeability", moment_tensor_2d_zz, pressure_units());
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
