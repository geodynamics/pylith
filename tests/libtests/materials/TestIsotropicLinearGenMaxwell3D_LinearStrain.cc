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

#include "TestIsotropicLinearGenMaxwell3D.hh" // Implementation of cases

#include "pylith/materials/IsotropicLinearGenMaxwell3D.hh" // USES IsotropicLinearGenMaxwell3D
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // :TEMPORARY: USES PYLITH_JOURNAL_ERROR

// forward declarations
namespace pylith {
    namespace materials {
        class TestIsotropicLinearGenMaxwell3D_LinearStrain;

        class TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP1;
        class TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP2;
        class TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP3;
        class TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP4;

        class TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ1;
        class TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ2;
        class TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ3;
        class TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ4;

    } // materials
} // pylith


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain :
    public pylith::materials::TestIsotropicLinearGenMaxwell3D {

    /// Spatial database user functions for auxiiliary subfields (includes derived fields).
    static const double SMALL;

    // Density
    static double density(const double x,
                          const double y,
                          const double z) {
        return 4000.0;
    } // density
	
    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units
	
    // Vs
    static double vs(const double x,
                     const double y,
                     const double z) {
        return 5600.0;
    } // vs
	
    static const char* vs_units(void) {
        return "m/s";
    } // vs_units
	
    // Vp
    static double vp(const double x,
                     const double y,
                     const double z) {
        return sqrt(3.0)*vs(x,y,z);
    } // vp

    static const char* vp_units(void) {
        return "m/s";
    } // vp_units

    // Viscosity_1
    static double viscosity_1(const double x,
                              const double y,
                              const double z) {
        return 7.91700159488e+19;
    } // viscosity_1

    // Viscosity_2
    static double viscosity_2(const double x,
                              const double y,
                              const double z) {
        return 7.91700159488e+18;
    } // viscosity_2

    // Viscosity_3
    static double viscosity_3(const double x,
                              const double y,
                              const double z) {
        return 7.91700159488e+17;
    } // viscosity_3

    static const char* viscosity_units(void) {
        return "Pa*s";
    } // viscosity__units

    // shear modulus
    static double shearModulus(const double x,
                               const double y,
                               const double z) {
        return density(x,y,z) * vs(x,y,z) * vs(x,y,z);
    } // shearModulus

    static const char* shearModulus_units(void) {
        return "Pa";
    } // shearModulus_units

    // bulk modulus
    static double bulkModulus(const double x,
                              const double y,
                              const double z) {
        return density(x,y,z)*(vp(x,y,z)*vp(x,y,z) - 4.0/3.0*vs(x,y,z)*vs(x,y,z));
    } // bulkModulus

    static const char* bulkModulus_units(void) {
        return "Pa";
    } // bulkModulus_units

    // Maxwell time_1
    static double maxwellTime_1(const double x,
                                const double y,
                                const double z) {
        return viscosity_1(x,y,z) / shearModulus(x,y,z);
    } // maxwellTime_1

    // Maxwell time_2
    static double maxwellTime_2(const double x,
                                const double y,
                                const double z) {
        return viscosity_2(x,y,z) / shearModulus(x,y,z);
    } // maxwellTime_2

    // Maxwell time_3
    static double maxwellTime_3(const double x,
                                const double y,
                                const double z) {
        return viscosity_3(x,y,z) / shearModulus(x,y,z);
    } // maxwellTime_3

    static const char* maxwellTime_units(void) {
        return "s";
    } // maxwellTime_units

    // shearModulusRatio_1
    static double shearModulusRatio_1(const double x,
                                      const double y,
                                      const double z) {
        return 0.4;
    } // shearModulusRatio_1

    // shearModulusRatio_2
    static double shearModulusRatio_2(const double x,
                                      const double y,
                                      const double z) {
        return 0.3;
    } // shearModulusRatio_2

    // shearModulusRatio_3
    static double shearModulusRatio_3(const double x,
                                      const double y,
                                      const double z) {
        return 0.2;
    } // shearModulusRatio_3

    struct AuxConstants {
        double a;
        double b;
        double c;
        double d;
        double e;
        double f;
        double g;
        double t;
        double dt;
    };
    static const AuxConstants constants;

    // Total strain

    static double totalStrain_xx(const double x,
                                 const double y,
                                 const double z) {
        return (2.0*constants.a*x + 2.0*constants.b*y + 2.0*constants.d*z) *
			(shearModulusRatio_1(x,y,z) * exp(-constants.t/maxwellTime_1(x,y,z)) +
			 shearModulusRatio_2(x,y,z) * exp(-constants.t/maxwellTime_2(x,y,z)) +
			 shearModulusRatio_3(x,y,z) * exp(-constants.t/maxwellTime_3(x,y,z)));
    } // totalStrain_xx

    static double totalStrain_yy(const double x,
                                 const double y,
                                 const double z) {
        return (2.0*constants.a*y + 2.0*constants.b*x + 2.0*constants.e*z) *
			(shearModulusRatio_1(x,y,z) * exp(-constants.t/maxwellTime_1(x,y,z)) +
			 shearModulusRatio_2(x,y,z) * exp(-constants.t/maxwellTime_2(x,y,z)) +
			 shearModulusRatio_3(x,y,z) * exp(-constants.t/maxwellTime_3(x,y,z)));
    } // totalStrain_yy

    static double totalStrain_zz(const double x,
                                 const double y,
                                 const double z) {
        return (2.0*constants.a*z + 2.0*constants.d*x + 2.0*constants.e*y) *
			(shearModulusRatio_1(x,y,z) * exp(-constants.t/maxwellTime_1(x,y,z)) +
			 shearModulusRatio_2(x,y,z) * exp(-constants.t/maxwellTime_2(x,y,z)) +
			 shearModulusRatio_3(x,y,z) * exp(-constants.t/maxwellTime_3(x,y,z)));
    } // totalStrain_zz

    static double totalStrain_xy(const double x,
                                 const double y,
                                 const double z) {
        return (constants.b * (x + y) + constants.c * (x + y) + z * (constants.d + constants.e)) *
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z))))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)));
    } // totalStrain_xy

    static double totalStrain_yz(const double x,
                                 const double y,
                                 const double z) {
        return (x*(constants.b + constants.d) + y*(constants.c + constants.e) + z*(constants.e + constants.f))*
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z))))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)));
    } // totalStrain_yz

    static double totalStrain_xz(const double x,
                                 const double y,
                                 const double z) {
        return (y*(constants.b + constants.e) + constants.d*(x + z) + constants.f*(x + z))*
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z))))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)));
    } // totalStrain_xz

    // Viscous strain 1

    static double viscousStrain_1_xx(const double x,
                                     const double y,
                                     const double z) {
        return 2.0*maxwellTime_1(x,y,z)*(exp(constants.t/maxwellTime_1(x,y,z)) - 1.0)*
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z))))*
			(constants.a*(2.0*x - y - z) + constants.b*(2.0*y - x) + constants.d*(2.0*z - x) -constants.e*(y + z))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 2.0*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/(3.0*constants.t);
    } // viscousStrain_1_xx

    static double viscousStrain_1_yy(const double x,
                                     const double y,
                                     const double z) {
        return -2.0*maxwellTime_1(x,y,z)*(exp(constants.t/maxwellTime_1(x,y,z)) - 1.0)*
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z))))*
			(constants.a*(x - 2.0*y + z) - constants.b*(2.0*x - y) + constants.d*(x + z) + constants.e*(y - 2.0*z))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 2.0*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/(3.0*constants.t);
    } // viscousStrain_1_yy

    static double viscousStrain_1_zz(const double x,
                                     const double y,
                                     const double z) {
        return -2.0*maxwellTime_1(x,y,z)*(exp(constants.t/maxwellTime_1(x,y,z)) - 1.0)*
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z))))*
			(constants.a*(x + y - 2.0*z) + constants.b*(x + y) -constants.d*(2.0*x - z) - constants.e*(2.0*y - z))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 2.0*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/(3.0*constants.t);
    } // viscousStrain_1_zz

    static double viscousStrain_1_xy(const double x,
                                     const double y,
                                     const double z) {
        return maxwellTime_1(x,y,z)*(exp(constants.t/maxwellTime_1(x,y,z)) - 1.0)*
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z))))*
			(constants.b*(x + y) + constants.c*(x + y) + z*(constants.d + constants.e))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 2.0*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/constants.t;
    } // viscousStrain_1_xy

    static double viscousStrain_1_yz(const double x,
                                     const double y,
                                     const double z) {
        return maxwellTime_1(x,y,z)*(exp(constants.t/maxwellTime_1(x,y,z)) - 1.0)*
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z))))*
			(x*(constants.b + constants.d) + y*(constants.c + constants.e) + z*(constants.e + constants.f))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 2.0*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/constants.t;
    } // viscousStrain_1_yz

    static double viscousStrain_1_xz(const double x,
                                     const double y,
                                     const double z) {
        return maxwellTime_1(x,y,z)*(exp(constants.t/maxwellTime_1(x,y,z)) - 1.0)*
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z))))*
			(y*(constants.b + constants.e) + constants.d*(x + z) + constants.f*(x + z))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 2.0*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/constants.t;
    } // viscousStrain_1_xz

    // Viscous strain 2

    static double viscousStrain_2_xx(const double x,
                                     const double y,
                                     const double z) {
        return 2.0*maxwellTime_2(x,y,z)*(exp(constants.t/maxwellTime_2(x,y,z)) - 1.0)*
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z))))*
			(constants.a*(2.0*x - y - z) + constants.b*(2.0*y - x) + constants.d*(2.0*z - x) -constants.e*(y + z))*
               exp(-constants.t*
                   (maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + 2.0*maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
					maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				   (maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/(3.0*constants.t);
    } // viscousStrain_2_xx

    static double viscousStrain_2_yy(const double x,
                                     const double y,
                                     const double z) {
        return -2.0*maxwellTime_2(x,y,z)*(exp(constants.t/maxwellTime_2(x,y,z)) - 1.0)*
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z))))*
			(constants.a*(x - 2.0*y + z) - constants.b*(2.0*x - y) + constants.d*(x + z) - constants.e*(2.0*z - y))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + 2.0*maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/(3.0*constants.t);
    } // viscousStrain_2_yy

    static double viscousStrain_2_zz(const double x,
                                     const double y,
                                     const double z) {
        return -2.0*maxwellTime_2(x,y,z)*(exp(constants.t/maxwellTime_2(x,y,z)) - 1.0) *
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)))) *
			(constants.a*(x + y - 2.0*z) + constants.b*(x + y) + constants.d*(z - 2.0*x) + constants.e*(z - 2.0*y))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + 2.0*maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/(3.0*constants.t);
    } // viscousStrain_2_zz

    static double viscousStrain_2_xy(const double x,
                                     const double y,
                                     const double z) {
        return maxwellTime_2(x,y,z)*(exp(constants.t/maxwellTime_2(x,y,z)) - 1.0) *
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)))) *
			(constants.b*(x + y) + constants.c*(x + y) + z*(constants.d + constants.e))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + 2.0*maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/constants.t;
    } // viscousStrain_2_xy

    static double viscousStrain_2_yz(const double x,
                                     const double y,
                                     const double z) {
        return maxwellTime_2(x,y,z)*(exp(constants.t/maxwellTime_2(x,y,z)) - 1.0) *
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)))) *
			(x*(constants.b + constants.d) + y*(constants.c + constants.e) + z*(constants.e + constants.f))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + 2.0*maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/constants.t;
    } // viscousStrain_2_yz

    static double viscousStrain_2_xz(const double x,
                                     const double y,
                                     const double z) {
        return maxwellTime_2(x,y,z)*(exp(constants.t/maxwellTime_2(x,y,z)) - 1.0) *
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)))) *
			(y*(constants.b + constants.e) + constants.d*(x + z) + constants.f*(x + z))*
			exp(-constants.t*
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + 2.0*maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/constants.t;
    } // viscousStrain_2_xz

    // Viscous strain 3

    static double viscousStrain_3_xx(const double x,
                                     const double y,
                                     const double z) {
        return 2.0*maxwellTime_3(x,y,z)*(exp(constants.t/maxwellTime_3(x,y,z)) - 1.0) *
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)))) *
			(constants.a*(2.0*x - y - z) + constants.b*(2.0*y - x) + constants.d*(2.0*z - x) -constants.e*(y + z))*
			exp(-constants.t*
				(2.0*maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/(3.0*constants.t);
    } // viscousStrain_3_xx

    static double viscousStrain_3_yy(const double x,
                                     const double y,
                                     const double z) {
        return -2.0*maxwellTime_3(x,y,z)*(exp(constants.t/maxwellTime_3(x,y,z)) - 1.0) *
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)))) *
			(constants.a*(x - 2.0*y + z) - constants.b*(2.0*x - y) + constants.d*(x + z) - constants.e*(2.0*z - y))*
			exp(-constants.t*
				(2.0*maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/(3.0*constants.t);
    } // viscousStrain_3_yy

    static double viscousStrain_3_zz(const double x,
                                     const double y,
                                     const double z) {
        return -2.0*maxwellTime_3(x,y,z)*(exp(constants.t/maxwellTime_3(x,y,z)) - 1.0) *
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)))) *
			(constants.a*(x + y - 2.0*z) + constants.b*(x + y) + constants.d*(z - 2.0*x) + constants.e*(z - 2.0*y))*
			exp(-constants.t*
				(2.0*maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/(3.0*constants.t);
    } // viscousStrain_3_zz

    static double viscousStrain_3_xy(const double x,
                                     const double y,
                                     const double z) {
        return maxwellTime_3(x,y,z)*(exp(constants.t/maxwellTime_3(x,y,z)) - 1.0) *
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)))) *
			(constants.b*(x + y) + constants.c*(x + y) + z*(constants.d + constants.e))*
			exp(-constants.t*
				(2.0*maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/constants.t;
    } // viscousStrain_3_xy

    static double viscousStrain_3_yz(const double x,
                                     const double y,
                                     const double z) {
        return maxwellTime_3(x,y,z)*(exp(constants.t/maxwellTime_3(x,y,z)) - 1.0) *
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)))) *
			(x*(constants.b + constants.d) + y*(constants.c + constants.e) + z*(constants.e + constants.f))*
			exp(-constants.t*
				(2.0*maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/constants.t;
    } // viscousStrain_3_yz

    static double viscousStrain_3_xz(const double x,
                                     const double y,
                                     const double z) {
        return maxwellTime_3(x,y,z)*(exp(constants.t/maxwellTime_3(x,y,z)) - 1.0) *
			(shearModulusRatio_1(x,y,z)*exp(constants.t*(maxwellTime_2(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_2(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_3(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z))) +
			 shearModulusRatio_3(x,y,z)*exp(constants.t*(maxwellTime_1(x,y,z) + maxwellTime_2(x,y,z))/
											(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)))) *
			(y*(constants.b + constants.e) + constants.d*(x + z) + constants.f*(x + z))*
			exp(-constants.t*
				(2.0*maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z) + maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z) +
				 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z))/
				(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)))/constants.t;
    } // viscousStrain_3_xz

    // Total strain for perturbed solution.

    static double totalStrainUpdate_xx(const double x,
                                       const double y,
                                       const double z) {
		return totalStrain_xx(x,y,z) + constants.g;
    } // totalStrainUpdate_xx

    static double totalStrainUpdate_yy(const double x,
                                       const double y,
                                       const double z) {
		return totalStrain_yy(x,y,z);
    } // totalStrainUpdate_yy

    static double totalStrainUpdate_zz(const double x,
                                       const double y,
                                       const double z) {
		return totalStrain_zz(x,y,z);
    } // totalStrainUpdate_zz

    static double totalStrainUpdate_xy(const double x,
                                       const double y,
                                       const double z) {
		return totalStrain_xy(x,y,z) + constants.g/2.0;
    } // totalStrainUpdate_xy

    static double totalStrainUpdate_yz(const double x,
                                       const double y,
                                       const double z) {
		return totalStrain_yz(x,y,z);
    } // totalStrainUpdate_yz

    static double totalStrainUpdate_xz(const double x,
                                       const double y,
                                       const double z) {
		return totalStrain_xz(x,y,z) + constants.g/2.0;
    } // totalStrainUpdate_xz

    // Values needed to compute viscous strain for perturbed solution.
    static double meanStrainT(const double x,
                              const double y,
                              const double z) {
        return (totalStrain_xx(x,y,z) + totalStrain_yy(x,y,z) + totalStrain_zz(x,y,z))/3.0;
    }

    static double devStrainT_xx(const double x,
                                const double y,
                                const double z) {
        return totalStrain_xx(x,y,z) - meanStrainT(x,y,z);
    }

    static double devStrainT_yy(const double x,
                                const double y,
                                const double z) {
        return totalStrain_yy(x,y,z) - meanStrainT(x,y,z);
    }

    static double devStrainT_zz(const double x,
                                const double y,
                                const double z) {
        return totalStrain_zz(x,y,z) - meanStrainT(x,y,z);
    }

    static double devStrainT_xy(const double x,
                                const double y,
                                const double z) {
        return totalStrain_xy(x,y,z);
    }

    static double devStrainT_yz(const double x,
                                const double y,
                                const double z) {
        return totalStrain_yz(x,y,z);
    }

    static double devStrainT_xz(const double x,
                                const double y,
                                const double z) {
        return totalStrain_xz(x,y,z);
    }

    static double meanStrainTplusDt(const double x,
                                    const double y,
                                    const double z) {
        return (totalStrainUpdate_xx(x,y,z) + totalStrainUpdate_yy(x,y,z) + totalStrainUpdate_zz(x,y,z))/3.0;
    }

    static double devStrainTplusDt_xx(const double x,
                                      const double y,
                                      const double z) {
        return totalStrainUpdate_xx(x,y,z) - meanStrainTplusDt(x,y,z);
    }

    static double devStrainTplusDt_yy(const double x,
                                      const double y,
                                      const double z) {
        return totalStrainUpdate_yy(x,y,z) - meanStrainTplusDt(x,y,z);
    }

    static double devStrainTplusDt_zz(const double x,
                                      const double y,
                                      const double z) {
        return totalStrainUpdate_zz(x,y,z) - meanStrainTplusDt(x,y,z);
    }

    static double devStrainTplusDt_xy(const double x,
                                      const double y,
                                      const double z) {
        return totalStrainUpdate_xy(x,y,z);
    }

    static double devStrainTplusDt_yz(const double x,
                                      const double y,
                                      const double z) {
        return totalStrainUpdate_yz(x,y,z);
    }

    static double devStrainTplusDt_xz(const double x,
                                      const double y,
                                      const double z) {
        return totalStrainUpdate_xz(x,y,z);
    }

    static double dq_1(const double x,
                       const double y,
                       const double z) {
        return maxwellTime_1(x,y,z) * (1.0 - exp(-constants.dt/maxwellTime_1(x,y,z)))/constants.dt;
    }

    static double dq_2(const double x,
                       const double y,
                       const double z) {
        return maxwellTime_2(x,y,z) * (1.0 - exp(-constants.dt/maxwellTime_2(x,y,z)))/constants.dt;
    }

    static double dq_3(const double x,
                       const double y,
                       const double z) {
        return maxwellTime_3(x,y,z) * (1.0 - exp(-constants.dt/maxwellTime_3(x,y,z)))/constants.dt;
    }

    // Viscous strain for perturbed solution.

    static double viscousStrainUpdate_1_xx(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_1(x,y,z))*viscousStrain_1_xx(x,y,z) +
			dq_1(x,y,z)*(devStrainTplusDt_xx(x,y,z) - devStrainT_xx(x,y,z));
    } // viscousStrainUpdate_1_xx

    static double viscousStrainUpdate_1_yy(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_1(x,y,z))*viscousStrain_1_yy(x,y,z) +
			dq_1(x,y,z)*(devStrainTplusDt_yy(x,y,z) - devStrainT_yy(x,y,z));
    } // viscousStrainUpdate_1_yy

    static double viscousStrainUpdate_1_zz(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_1(x,y,z))*viscousStrain_1_zz(x,y,z) +
			dq_1(x,y,z)*(devStrainTplusDt_zz(x,y,z) - devStrainT_zz(x,y,z));
    } // viscousStrainUpdate_1_zz

    static double viscousStrainUpdate_1_xy(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_1(x,y,z))*viscousStrain_1_xy(x,y,z) +
			dq_1(x,y,z)*(devStrainTplusDt_xy(x,y,z) - devStrainT_xy(x,y,z));
    } // viscousStrainUpdate_1_xy

    static double viscousStrainUpdate_1_yz(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_1(x,y,z))*viscousStrain_1_yz(x,y,z) +
			dq_1(x,y,z)*(devStrainTplusDt_yz(x,y,z) - devStrainT_yz(x,y,z));
    } // viscousStrainUpdate_1_yz

    static double viscousStrainUpdate_1_xz(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_1(x,y,z))*viscousStrain_1_xz(x,y,z) +
			dq_1(x,y,z)*(devStrainTplusDt_xz(x,y,z) - devStrainT_xz(x,y,z));
    } // viscousStrainUpdate_1_xz

    static double viscousStrainUpdate_2_xx(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_2(x,y,z))*viscousStrain_2_xx(x,y,z) +
			dq_2(x,y,z)*(devStrainTplusDt_xx(x,y,z) - devStrainT_xx(x,y,z));
    } // viscousStrainUpdate_2_xx

    static double viscousStrainUpdate_2_yy(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_2(x,y,z))*viscousStrain_2_yy(x,y,z) +
			dq_2(x,y,z)*(devStrainTplusDt_yy(x,y,z) - devStrainT_yy(x,y,z));
    } // viscousStrainUpdate_2_yy

    static double viscousStrainUpdate_2_zz(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_2(x,y,z))*viscousStrain_2_zz(x,y,z) +
			dq_2(x,y,z)*(devStrainTplusDt_zz(x,y,z) - devStrainT_zz(x,y,z));
    } // viscousStrainUpdate_2_zz

    static double viscousStrainUpdate_2_xy(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_2(x,y,z))*viscousStrain_2_xy(x,y,z) +
			dq_2(x,y,z)*(devStrainTplusDt_xy(x,y,z) - devStrainT_xy(x,y,z));
    } // viscousStrainUpdate_2_xy

    static double viscousStrainUpdate_2_yz(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_2(x,y,z))*viscousStrain_2_yz(x,y,z) +
			dq_2(x,y,z)*(devStrainTplusDt_yz(x,y,z) - devStrainT_yz(x,y,z));
    } // viscousStrainUpdate_2_yz

    static double viscousStrainUpdate_2_xz(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_2(x,y,z))*viscousStrain_2_xz(x,y,z) +
			dq_2(x,y,z)*(devStrainTplusDt_xz(x,y,z) - devStrainT_xz(x,y,z));
    } // viscousStrainUpdate_2_xz

    static double viscousStrainUpdate_3_xx(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_3(x,y,z))*viscousStrain_3_xx(x,y,z) +
			dq_3(x,y,z)*(devStrainTplusDt_xx(x,y,z) - devStrainT_xx(x,y,z));
    } // viscousStrainUpdate_3_xx

    static double viscousStrainUpdate_3_yy(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_3(x,y,z))*viscousStrain_3_yy(x,y,z) +
			dq_3(x,y,z)*(devStrainTplusDt_yy(x,y,z) - devStrainT_yy(x,y,z));
    } // viscousStrainUpdate_3_yy

    static double viscousStrainUpdate_3_zz(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_3(x,y,z))*viscousStrain_3_zz(x,y,z) +
			dq_3(x,y,z)*(devStrainTplusDt_zz(x,y,z) - devStrainT_zz(x,y,z));
    } // viscousStrainUpdate_3_zz

    static double viscousStrainUpdate_3_xy(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_3(x,y,z))*viscousStrain_3_xy(x,y,z) +
			dq_3(x,y,z)*(devStrainTplusDt_xy(x,y,z) - devStrainT_xy(x,y,z));
    } // viscousStrainUpdate_3_xy

    static double viscousStrainUpdate_3_yz(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_3(x,y,z))*viscousStrain_3_yz(x,y,z) +
			dq_3(x,y,z)*(devStrainTplusDt_yz(x,y,z) - devStrainT_yz(x,y,z));
    } // viscousStrainUpdate_3_yz

    static double viscousStrainUpdate_3_xz(const double x,
                                           const double y,
                                           const double z) {
        return exp(-constants.dt/maxwellTime_3(x,y,z))*viscousStrain_3_xz(x,y,z) +
			dq_3(x,y,z)*(devStrainTplusDt_xz(x,y,z) - devStrainT_xz(x,y,z));
    } // viscousStrainUpdate_3_xz

    // Body force
    static double bodyforce_x(const double x,
                              const double y,
                              const double z) {
        return 6.0*(constants.a + constants.b + constants.d)*
			bulkModulus(x,y,z)*shearModulusRatio_1(x,y,z)*exp(-constants.t/maxwellTime_1(x,y,z)) +
			6.0*(constants.a + constants.b + constants.d)*
			bulkModulus(x,y,z)*shearModulusRatio_2(x,y,z)*exp(-constants.t/maxwellTime_2(x,y,z)) +
			6.0*(constants.a + constants.b + constants.d)*
			bulkModulus(x,y,z)*shearModulusRatio_3(x,y,z)*exp(-constants.t/maxwellTime_3(x,y,z));

    } // bodyforce_x

    static double bodyforce_y(const double x,
                              const double y,
                              const double z) {
        return 6.0*(constants.a + constants.b + constants.e)*
			bulkModulus(x,y,z)*shearModulusRatio_1(x,y,z)*exp(-constants.t/maxwellTime_1(x,y,z)) +
			6.0*(constants.a + constants.b + constants.e)*
			bulkModulus(x,y,z)*shearModulusRatio_2(x,y,z)*exp(-constants.t/maxwellTime_2(x,y,z)) +
			6.0*(constants.a + constants.b + constants.e)*
			bulkModulus(x,y,z)*shearModulusRatio_3(x,y,z)*exp(-constants.t/maxwellTime_3(x,y,z));

    } // bodyforce_y

    static double bodyforce_z(const double x,
                              const double y,
                              const double z) {
        return 6.0*(constants.a + constants.d + constants.e)*
			bulkModulus(x,y,z)*shearModulusRatio_1(x,y,z)*exp(-constants.t/maxwellTime_1(x,y,z)) +
			6.0*(constants.a + constants.d + constants.e)*
			bulkModulus(x,y,z)*shearModulusRatio_2(x,y,z)*exp(-constants.t/maxwellTime_2(x,y,z)) +
			6.0*(constants.a + constants.d + constants.e)*
			bulkModulus(x,y,z)*shearModulusRatio_3(x,y,z)*exp(-constants.t/maxwellTime_3(x,y,z));

    } // bodyforce_z

    static const char* bodyforce_units(void) {
        return "kg*m/s**2";
    }

    // Spatial database user functions for solution subfields.

    // Displacement
    static double disp_x(const double x,
                         const double y,
                         const double z) {
        return (shearModulusRatio_1(x,y,z)*exp(constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_2(x,y,z))) +
                shearModulusRatio_2(x,y,z)*exp(constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_1(x,y,z))) +
                shearModulusRatio_3(x,y,z)*exp(constants.t*(1.0/maxwellTime_2(x,y,z) + 1.0/maxwellTime_1(x,y,z))))*
			(constants.a*x*x + 2.0*constants.b*x*y + constants.c*y*y + 2.0*constants.d*x*z + 2.0*constants.e*y*z +
			 constants.f*z*z)*
			exp(-constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_2(x,y,z) + 1.0/maxwellTime_1(x,y,z)));

    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double z) {
        return (shearModulusRatio_1(x,y,z)*exp(constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_2(x,y,z))) +
                shearModulusRatio_2(x,y,z)*exp(constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_1(x,y,z))) +
                shearModulusRatio_3(x,y,z)*exp(constants.t*(1.0/maxwellTime_2(x,y,z) + 1.0/maxwellTime_1(x,y,z))))*
			(constants.a*y*y + 2.0*constants.b*x*y + constants.c*x*x + 2.0*constants.d*x*z + 2.0*constants.e*y*z +
			 constants.f*z*z)*
			exp(-constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_2(x,y,z) + 1.0/maxwellTime_1(x,y,z)));

    } // disp_y

    static double disp_z(const double x,
                         const double y,
                         const double z) {
        return (shearModulusRatio_1(x,y,z)*exp(constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_2(x,y,z))) +
                shearModulusRatio_2(x,y,z)*exp(constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_1(x,y,z))) +
                shearModulusRatio_3(x,y,z)*exp(constants.t*(1.0/maxwellTime_2(x,y,z) + 1.0/maxwellTime_1(x,y,z))))*
			(constants.a*z*z + 2.0*constants.b*x*y + constants.c*y*y + 2.0*constants.d*x*z + 2.0*constants.e*y*z +
			 constants.f*x*x)*
			exp(-constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_2(x,y,z) + 1.0/maxwellTime_1(x,y,z)));

    } // disp_z

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    static double disp_dot_x(const double x,
                             const double y,
                             const double z) {
        return -(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*shearModulusRatio_3(x,y,z)*
                 exp(constants.t*(1.0/maxwellTime_2(x,y,z) + 1.0/maxwellTime_1(x,y,z))) +
                 maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z)*shearModulusRatio_2(x,y,z)*
                 exp(constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_1(x,y,z))) +
                 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)*shearModulusRatio_1(x,y,z)*
                 exp(constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_2(x,y,z))))*
			(constants.a*x*x + 2.0*constants.b*x*y + constants.c*y*y + 2.0*constants.d*x*z + 2.0*constants.e*y*z +
			 constants.f*z*z)*
			exp(-constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_2(x,y,z) + 1.0/maxwellTime_1(x,y,z)))/
			(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z));

    } // disp_dot_x

    static double disp_dot_y(const double x,
                             const double y,
                             const double z) {
        return -(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*shearModulusRatio_3(x,y,z)*
                 exp(constants.t*(1.0/maxwellTime_2(x,y,z) + 1.0/maxwellTime_1(x,y,z))) +
                 maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z)*shearModulusRatio_2(x,y,z)*
                 exp(constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_1(x,y,z))) +
                 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)*shearModulusRatio_1(x,y,z)*
                 exp(constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_2(x,y,z))))*
			(constants.a*y*y + 2.0*constants.b*x*y + constants.c*x*x + 2.0*constants.d*x*z + 2.0*constants.e*y*z +
			 constants.f*z*z)*
			exp(-constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_2(x,y,z) + 1.0/maxwellTime_1(x,y,z)))/
			(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z));

    } // disp_dot_y

    static double disp_dot_z(const double x,
                             const double y,
                             const double z) {
        return -(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*shearModulusRatio_3(x,y,z)*
                 exp(constants.t*(1.0/maxwellTime_2(x,y,z) + 1.0/maxwellTime_1(x,y,z))) +
                 maxwellTime_1(x,y,z)*maxwellTime_3(x,y,z)*shearModulusRatio_2(x,y,z)*
                 exp(constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_1(x,y,z))) +
                 maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z)*shearModulusRatio_1(x,y,z)*
                 exp(constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_2(x,y,z))))*
			(constants.a*z*z + 2.0*constants.b*x*y + constants.c*y*y + 2.0*constants.d*x*z + 2.0*constants.e*y*z +
			 constants.f*x*x)*
			exp(-constants.t*(1.0/maxwellTime_3(x,y,z) + 1.0/maxwellTime_2(x,y,z) + 1.0/maxwellTime_1(x,y,z)))/
			(maxwellTime_1(x,y,z)*maxwellTime_2(x,y,z)*maxwellTime_3(x,y,z));

    } // disp_dot_z

    static const char* disp_dot_units(void) {
        return "m/s";
    } // disp_dot_units

    // Displacement + perturbation
    static double disp_perturb_x(const double x,
                                 const double y,
                                 const double z) {
        return disp_x(x,y,z) + constants.g * x;
    } // disp_perturb_x

    static double disp_perturb_y(const double x,
                                 const double y,
                                 const double z) {
        return disp_y(x,y,z) + constants.g * x;
    } // disp_perturb_y

    static double disp_perturb_z(const double x,
                                 const double y,
                                 const double z) {
        return disp_z(x,y,z) + constants.g * x;
    } // disp_perturb_z

protected:
    void setUp(void) {
        TestIsotropicLinearGenMaxwell3D::setUp();
        _mydata = new TestIsotropicLinearGenMaxwell3D_Data(); CPPUNIT_ASSERT(_mydata);

        // dimension set in base class.
        // meshFilename set in derived class.
        _mydata->boundaryLabel = "boundary";

        CPPUNIT_ASSERT(_mydata->normalizer);
        _mydata->normalizer->setLengthScale(1.0e+03);
        _mydata->normalizer->setTimeScale(2.0e+7);
        _mydata->normalizer->setDensityScale(3.0e+3);
        _mydata->normalizer->setPressureScale(1.25e+11);

        _mydata->t = constants.t/_mydata->normalizer->getTimeScale();
        _mydata->dt = constants.dt/_mydata->normalizer->getTimeScale();
        _mydata->s_tshift = 1.0 / _mydata->dt;

        // solnDiscretizations set in derived class.

        _mydata->numAuxSubfields = 8;
        static const char* _auxSubfields[8] =
        {"density", "shear_modulus", "bulk_modulus",
         "maxwell_time", "shear_modulus_ratio", "viscous_strain",
         "total_strain", "body_force"};
        _mydata->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // maxwell_time
            pylith::topology::Field::Discretization(0, 1), // shear_modulus_ratio
            pylith::topology::Field::Discretization(1, 1, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 1, false), // total_strain
            pylith::topology::Field::Discretization(0, 1), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        CPPUNIT_ASSERT(_mydata->auxDB);
        _mydata->auxDB->addValue("density", density, density_units());
        _mydata->auxDB->addValue("vp", vp, vp_units());
        _mydata->auxDB->addValue("vs", vs, vs_units());
        _mydata->auxDB->addValue("shear_modulus", shearModulus, shearModulus_units());
        _mydata->auxDB->addValue("bulk_modulus", bulkModulus, bulkModulus_units());
        _mydata->auxDB->addValue("viscosity_1", viscosity_1, viscosity_units());
        _mydata->auxDB->addValue("viscosity_2", viscosity_2, viscosity_units());
        _mydata->auxDB->addValue("viscosity_3", viscosity_3, viscosity_units());
        _mydata->auxDB->addValue("maxwell_time_1", maxwellTime_1, maxwellTime_units());
        _mydata->auxDB->addValue("maxwell_time_2", maxwellTime_2, maxwellTime_units());
        _mydata->auxDB->addValue("maxwell_time_3", maxwellTime_3, maxwellTime_units());
        _mydata->auxDB->addValue("shear_modulus_ratio_1", shearModulusRatio_1, "none");
        _mydata->auxDB->addValue("shear_modulus_ratio_2", shearModulusRatio_2, "none");
        _mydata->auxDB->addValue("shear_modulus_ratio_3", shearModulusRatio_3, "none");
        _mydata->auxDB->addValue("viscous_strain_1_xx", viscousStrain_1_xx, "none");
        _mydata->auxDB->addValue("viscous_strain_1_yy", viscousStrain_1_yy, "none");
        _mydata->auxDB->addValue("viscous_strain_1_zz", viscousStrain_1_zz, "none");
        _mydata->auxDB->addValue("viscous_strain_1_xy", viscousStrain_1_xy, "none");
        _mydata->auxDB->addValue("viscous_strain_1_yz", viscousStrain_1_yz, "none");
        _mydata->auxDB->addValue("viscous_strain_1_xz", viscousStrain_1_xz, "none");
        _mydata->auxDB->addValue("viscous_strain_2_xx", viscousStrain_2_xx, "none");
        _mydata->auxDB->addValue("viscous_strain_2_yy", viscousStrain_2_yy, "none");
        _mydata->auxDB->addValue("viscous_strain_2_zz", viscousStrain_2_zz, "none");
        _mydata->auxDB->addValue("viscous_strain_2_xy", viscousStrain_2_xy, "none");
        _mydata->auxDB->addValue("viscous_strain_2_yz", viscousStrain_2_yz, "none");
        _mydata->auxDB->addValue("viscous_strain_2_xz", viscousStrain_2_xz, "none");
        _mydata->auxDB->addValue("viscous_strain_3_xx", viscousStrain_3_xx, "none");
        _mydata->auxDB->addValue("viscous_strain_3_yy", viscousStrain_3_yy, "none");
        _mydata->auxDB->addValue("viscous_strain_3_zz", viscousStrain_3_zz, "none");
        _mydata->auxDB->addValue("viscous_strain_3_xy", viscousStrain_3_xy, "none");
        _mydata->auxDB->addValue("viscous_strain_3_yz", viscousStrain_3_yz, "none");
        _mydata->auxDB->addValue("viscous_strain_3_xz", viscousStrain_3_xz, "none");
        _mydata->auxDB->addValue("total_strain_xx", totalStrain_xx, "none");
        _mydata->auxDB->addValue("total_strain_yy", totalStrain_yy, "none");
        _mydata->auxDB->addValue("total_strain_zz", totalStrain_zz, "none");
        _mydata->auxDB->addValue("total_strain_xy", totalStrain_xy, "none");
        _mydata->auxDB->addValue("total_strain_yz", totalStrain_yz, "none");
        _mydata->auxDB->addValue("total_strain_xz", totalStrain_xz, "none");
        _mydata->auxDB->addValue("body_force_x", bodyforce_x, bodyforce_units());
        _mydata->auxDB->addValue("body_force_y", bodyforce_y, bodyforce_units());
        _mydata->auxDB->addValue("body_force_z", bodyforce_z, bodyforce_units());

        CPPUNIT_ASSERT(_mydata->auxUpdateDB);
        _mydata->auxUpdateDB->addValue("density", density, density_units());
        _mydata->auxUpdateDB->addValue("vp", vp, vp_units());
        _mydata->auxUpdateDB->addValue("vs", vs, vs_units());
        _mydata->auxUpdateDB->addValue("shear_modulus", shearModulus, shearModulus_units());
        _mydata->auxUpdateDB->addValue("bulk_modulus", bulkModulus, bulkModulus_units());
        _mydata->auxUpdateDB->addValue("viscosity_1", viscosity_1, viscosity_units());
        _mydata->auxUpdateDB->addValue("viscosity_2", viscosity_2, viscosity_units());
        _mydata->auxUpdateDB->addValue("viscosity_3", viscosity_3, viscosity_units());
        _mydata->auxUpdateDB->addValue("maxwell_time_1", maxwellTime_1, maxwellTime_units());
        _mydata->auxUpdateDB->addValue("maxwell_time_2", maxwellTime_2, maxwellTime_units());
        _mydata->auxUpdateDB->addValue("maxwell_time_3", maxwellTime_3, maxwellTime_units());
        _mydata->auxUpdateDB->addValue("shear_modulus_ratio_1", shearModulusRatio_1, "none");
        _mydata->auxUpdateDB->addValue("shear_modulus_ratio_2", shearModulusRatio_2, "none");
        _mydata->auxUpdateDB->addValue("shear_modulus_ratio_3", shearModulusRatio_3, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_1_xx", viscousStrainUpdate_1_xx, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_1_yy", viscousStrainUpdate_1_yy, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_1_zz", viscousStrainUpdate_1_zz, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_1_xy", viscousStrainUpdate_1_xy, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_1_yz", viscousStrainUpdate_1_yz, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_1_xz", viscousStrainUpdate_1_xz, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_2_xx", viscousStrainUpdate_2_xx, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_2_yy", viscousStrainUpdate_2_yy, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_2_zz", viscousStrainUpdate_2_zz, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_2_xy", viscousStrainUpdate_2_xy, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_2_yz", viscousStrainUpdate_2_yz, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_2_xz", viscousStrainUpdate_2_xz, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_3_xx", viscousStrainUpdate_3_xx, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_3_yy", viscousStrainUpdate_3_yy, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_3_zz", viscousStrainUpdate_3_zz, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_3_xy", viscousStrainUpdate_3_xy, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_3_yz", viscousStrainUpdate_3_yz, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_3_xz", viscousStrainUpdate_3_xz, "none");
        _mydata->auxUpdateDB->addValue("total_strain_xx", totalStrainUpdate_xx, "none");
        _mydata->auxUpdateDB->addValue("total_strain_yy", totalStrainUpdate_yy, "none");
        _mydata->auxUpdateDB->addValue("total_strain_zz", totalStrainUpdate_zz, "none");
        _mydata->auxUpdateDB->addValue("total_strain_xy", totalStrainUpdate_xy, "none");
        _mydata->auxUpdateDB->addValue("total_strain_yz", totalStrainUpdate_yz, "none");
        _mydata->auxUpdateDB->addValue("total_strain_xz", totalStrainUpdate_xz, "none");
        _mydata->auxUpdateDB->addValue("body_force_x", bodyforce_x, bodyforce_units());
        _mydata->auxUpdateDB->addValue("body_force_y", bodyforce_y, bodyforce_units());
        _mydata->auxUpdateDB->addValue("body_force_z", bodyforce_z, bodyforce_units());

        CPPUNIT_ASSERT(_mydata->solnDB);
        _mydata->solnDB->addValue("displacement_x", disp_x, disp_units());
        _mydata->solnDB->addValue("displacement_y", disp_y, disp_units());
        _mydata->solnDB->addValue("displacement_z", disp_z, disp_units());
        _mydata->solnDB->addValue("displacement_dot_x", disp_dot_x, disp_dot_units());
        _mydata->solnDB->addValue("displacement_dot_y", disp_dot_y, disp_dot_units());
        _mydata->solnDB->addValue("displacement_dot_z", disp_dot_z, disp_dot_units());

        CPPUNIT_ASSERT(_mydata->perturbDB);
        _mydata->perturbDB->addValue("displacement_x", disp_perturb_x, disp_units());
        _mydata->perturbDB->addValue("displacement_y", disp_perturb_y, disp_units());
        _mydata->perturbDB->addValue("displacement_z", disp_perturb_z, disp_units());
        _mydata->perturbDB->addValue("displacement_dot_x", disp_dot_x, disp_dot_units());
        _mydata->perturbDB->addValue("displacement_dot_y", disp_dot_y, disp_dot_units());
        _mydata->perturbDB->addValue("displacement_dot_z", disp_dot_z, disp_dot_units());

        CPPUNIT_ASSERT(_mymaterial);
        _mymaterial->useInertia(false);
        _mymaterial->useBodyForce(true);
        _mymaterial->useReferenceState(false);

        _mymaterial->setLabel("Isotropic Linear Generalized Maxwell 3D");
        _mymaterial->id(24);

    } // setUp

}; // TestIsotropicLinearGenMaxwell3D_LinearStrain
const double pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain::SMALL = 1.0e-5;

const pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain::AuxConstants pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain::constants = {
    1.0e-7, // a
    2.5e-7, // b
    3.0e-7, // c
    3.5e-7, // d
    4.0e-7, // e
    4.5e-7, // f
    9.0e-8, // g
    5.0e+7, // t
    5.0e+7, // dt
};


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP1 :
    public pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP1,
                           TestIsotropicLinearGenMaxwell3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearGenMaxwell3D_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tet_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP1
// Do not use this test because TetP1 cannot represent a linear strain field.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP1);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP2 :
    public pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP2,
                           TestIsotropicLinearGenMaxwell3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearGenMaxwell3D_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tet_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // maxwell_time
            pylith::topology::Field::Discretization(0, 2), // shear_modulus_ratio
            pylith::topology::Field::Discretization(1, 2, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 2, false), // total_strain
            pylith::topology::Field::Discretization(0, 2), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP2);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP3 :
    public pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP3,
                           TestIsotropicLinearGenMaxwell3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearGenMaxwell3D_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tet_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // maxwell_time
            pylith::topology::Field::Discretization(0, 3), // shear_modulus_ratio
            pylith::topology::Field::Discretization(1, 3, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 3, false), // total_strain
            pylith::topology::Field::Discretization(0, 3), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP3
// Remove this test for now until higher order integration is done properly.
//CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP3);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP4 :
    public pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP4,
                           TestIsotropicLinearGenMaxwell3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearGenMaxwell3D_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tet_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // maxwell_time
            pylith::topology::Field::Discretization(0, 4), // shear_modulus_ratio
            pylith::topology::Field::Discretization(1, 4, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 4, false), // total_strain
            pylith::topology::Field::Discretization(0, 4), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP4
// Remove this test for now until higher order integration is done properly.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_TetP4);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ1 :
    public pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ1,
                           TestIsotropicLinearGenMaxwell3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearGenMaxwell3D_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/hex_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // maxwell_time
            pylith::topology::Field::Discretization(0, 1), // shear_modulus_ratio
            pylith::topology::Field::Discretization(1, 1, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 1, false), // total_strain
            pylith::topology::Field::Discretization(0, 1), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ1
// Do not use this test because TetP1 cannot represent a linear strain field.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ1);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ2 :
    public pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ2,
                           TestIsotropicLinearGenMaxwell3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearGenMaxwell3D_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/hex_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // maxwell_time
            pylith::topology::Field::Discretization(0, 2), // shear_modulus_ratio
            pylith::topology::Field::Discretization(1, 2, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 2, false), // total_strain
            pylith::topology::Field::Discretization(0, 2), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ2);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ3 :
    public pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ3,
                           TestIsotropicLinearGenMaxwell3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearGenMaxwell3D_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/hex_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // maxwell_time
            pylith::topology::Field::Discretization(0, 3), // shear_modulus_ratio
            pylith::topology::Field::Discretization(1, 3, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 3, false), // total_strain
            pylith::topology::Field::Discretization(0, 3), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ3
// Remove this test for now until higher order integration is done properly.
//CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ3);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ4 :
    public pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ4,
                           TestIsotropicLinearGenMaxwell3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearGenMaxwell3D_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/hex_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // maxwell_time
            pylith::topology::Field::Discretization(0, 4), // shear_modulus_ratio
            pylith::topology::Field::Discretization(1, 4, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 4, false), // total_strain
            pylith::topology::Field::Discretization(0, 4), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ4
// Remove this test for now until higher order integration is done properly.
//CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwell3D_LinearStrain_HexQ4);


// End of file
