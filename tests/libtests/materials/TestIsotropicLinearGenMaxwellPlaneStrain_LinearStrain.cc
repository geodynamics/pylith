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

#include "TestIsotropicLinearGenMaxwellPlaneStrain.hh" // Implementation of cases

#include "pylith/materials/IsotropicLinearGenMaxwellPlaneStrain.hh" // USES IsotropicLinearGenMaxwellPlaneStrain
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // :TEMPORARY: USES PYLITH_JOURNAL_ERROR

// forward declarations
namespace pylith {
    namespace materials {
        class TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain;

        class TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP1;
        class TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP2;
        class TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP3;
        class TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP4;

        class TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ1;
        class TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ2;
        class TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ3;
        class TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ4;

    } // materials
} // pylith


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain :
	public pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain {

	/// Spatial database user functions for auxiiliary subfields (includes derived fields).
	static const double SMALL;

	// Density
	static double density(const double x,
						  const double y) {
		return 4000.0;
	} // density
	static const char* density_units(void) {
		return "kg/m**3";
	} // density_units

	// Vs
	static double vs(const double x,
					 const double y) {
		return 5600.0;
	} // vs
	static const char* vs_units(void) {
		return "m/s";
	} // vs_units

	// Vp
	static double vp(const double x,
					 const double y) {
		return sqrt(3.0)*vs(x,y);
	} // vp
	static const char* vp_units(void) {
		return "m/s";
	} // vp_units

	// Viscosity_1
	static double viscosity_1(const double x,
							  const double y) {
		return 7.91700159488e+19;
	} // viscosity_1

	// Viscosity_2
	static double viscosity_2(const double x,
							  const double y) {
		return 7.91700159488e+18;
	} // viscosity_2
	
	// Viscosity_3
	static double viscosity_3(const double x,
							  const double y) {
		return 7.91700159488e+17;
	} // viscosity_3
	
	static const char* viscosity_units(void) {
		return "Pa*s";
	} // viscosity__units
	
	// shear modulus
	static double shearModulus(const double x,
							   const double y) {
		return density(x,y) * vs(x,y) * vs(x,y);
	} // shearModulus
	static const char* shearModulus_units(void) {
		return "Pa";
	} // shearModulus_units
	
	// bulk modulus
	static double bulkModulus(const double x,
							  const double y) {
		return density(x,y)*(vp(x,y)*vp(x,y) - 4.0/3.0*vs(x,y)*vs(x,y));
	} // bulkModulus
	static const char* bulkModulus_units(void) {
		return "Pa";
	} // bulkModulus_units
	
	// Maxwell time_1
	static double maxwellTime_1(const double x,
								const double y) {
		return viscosity_1(x,y) / shearModulus(x,y);
	} // maxwellTime_1
	
	// Maxwell time_2
	static double maxwellTime_2(const double x,
								const double y) {
		return viscosity_2(x,y) / shearModulus(x,y);
	} // maxwellTime_2
	
	// Maxwell time_3
	static double maxwellTime_3(const double x,
								const double y) {
		return viscosity_3(x,y) / shearModulus(x,y);
	} // maxwellTime_3
	
	static const char* maxwellTime_units(void) {
		return "s";
	} // maxwellTime_units
	
	// shearModulusRatio_1
	static double shearModulusRatio_1(const double x,
									  const double y) {
		return 0.4;
	} // shearModulusRatio_1
	
	// shearModulusRatio_2
	static double shearModulusRatio_2(const double x,
									  const double y) {
		return 0.3;
	} // shearModulusRatio_2
	
	// shearModulusRatio_3
	static double shearModulusRatio_3(const double x,
									  const double y) {
		return 0.2;
	} // shearModulusRatio_3
	
	struct AuxConstants {
		double a;
		double b;
		double c;
		double d;
		double t;
		double dt;
	};
	static const AuxConstants constants;
	
	// Total strain
	
	static double totalStrain_xx(const double x,
								 const double y) {
		return (2.0*constants.a*x + 2.0*constants.b*y) *
			(shearModulusRatio_1(x,y) * exp(-constants.t/maxwellTime_1(x,y)) +
			 shearModulusRatio_2(x,y) * exp(-constants.t/maxwellTime_2(x,y)) +
			 shearModulusRatio_3(x,y) * exp(-constants.t/maxwellTime_3(x,y)));
	} // totalStrain_xx
	
	static double totalStrain_yy(const double x,
								 const double y) {
		return (2.0*constants.a*y + 2.0*constants.b*x) *
			(shearModulusRatio_1(x,y) * exp(-constants.t/maxwellTime_1(x,y)) +
			 shearModulusRatio_2(x,y) * exp(-constants.t/maxwellTime_2(x,y)) +
			 shearModulusRatio_3(x,y) * exp(-constants.t/maxwellTime_3(x,y)));
	} // totalStrain_yy
	
	static double totalStrain_zz(const double x,
								 const double y) {
		return 0.0;
	} // totalStrain_zz
	
	static double totalStrain_xy(const double x,
								 const double y) {
		return (shearModulusRatio_1(x,y)*exp(constants.t*(maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y)*maxwellTime_3(x,y))) +
				shearModulusRatio_2(x,y)*exp(constants.t*(maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y)*maxwellTime_3(x,y))) +
				shearModulusRatio_3(x,y)*exp(constants.t*(maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y)*maxwellTime_2(x,y))))*
			(constants.b * (x + y) + constants.c * (x + y)) *
			exp(-constants.t*(maxwellTime_1(x,y)*maxwellTime_2(x,y) + maxwellTime_1(x,y)*maxwellTime_3(x,y) + maxwellTime_2(x,y)*maxwellTime_3(x,y))/(maxwellTime_1(x,y)*maxwellTime_2(x,y)*maxwellTime_3(x,y)));
	} // totalStrain_xy
	
	// Viscous strain 1
	
	static double viscousStrain_1_xx(const double x,
									 const double y) {
		return 2.0 * maxwellTime_1(x,y) * (exp(constants.t/maxwellTime_1(x,y)) - 1) *
			(shearModulusRatio_1(x,y) * exp(constants.t * (maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_2(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_3(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y)))) *
			(constants.a * (2.0 * x - y) + constants.b * (2.0 * y - x)) *
			exp(-constants.t * (maxwellTime_1(x,y) * maxwellTime_2(x,y) + maxwellTime_1(x,y) * maxwellTime_3(x,y) + 2.0 * maxwellTime_2(x,y) * maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y) * maxwellTime_3(x,y)))/(3.0 * constants.t);
	} // viscousStrain_1_xx
	
	static double viscousStrain_1_yy(const double x,
									 const double y) {
		return -2.0 * maxwellTime_1(x,y) * (exp(constants.t/maxwellTime_1(x,y)) - 1.0) *
			(shearModulusRatio_1(x,y) * exp(constants.t * (maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_2(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_3(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y)))) *
			(constants.a * (x - 2.0 *y) - constants.b * (2.0 * x - y)) *
			exp(-constants.t * (maxwellTime_1(x,y) * maxwellTime_2(x,y) + maxwellTime_1(x,y) * maxwellTime_3(x,y) + 2.0 * maxwellTime_2(x,y) * maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y) * maxwellTime_3(x,y)))/(3.0 * constants.t);
	} // viscousStrain_1_yy
	
	static double viscousStrain_1_zz(const double x,
									 const double y) {
		return -2.0 * maxwellTime_1(x,y) * (exp(constants.t/maxwellTime_1(x,y)) - 1.0) *
			(shearModulusRatio_1(x,y) * exp(constants.t * (maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_2(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_3(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y)))) *
			(constants.a * (x + y) + constants.b * (x + y)) *
			exp(-constants.t * (maxwellTime_1(x,y) * maxwellTime_2(x,y) + maxwellTime_1(x,y) * maxwellTime_3(x,y) + 2.0 * maxwellTime_2(x,y) * maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y) * maxwellTime_3(x,y)))/(3.0 * constants.t);
	} // viscousStrain_1_zz

	static double viscousStrain_1_xy(const double x,
									 const double y) {
		return maxwellTime_1(x,y) * (exp(constants.t/maxwellTime_1(x,y)) - 1.0) *
			(shearModulusRatio_1(x,y) * exp(constants.t * (maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_2(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_3(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y)))) *
			(constants.b * (x + y) + constants.c * (x + y)) *
			exp(-constants.t * (maxwellTime_1(x,y) * maxwellTime_2(x,y) + maxwellTime_1(x,y) * maxwellTime_3(x,y) + 2.0 * maxwellTime_2(x,y) * maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y) * maxwellTime_3(x,y)))/constants.t;
	} // viscousStrain_1_xy
	
	// Viscous strain 2
	
	static double viscousStrain_2_xx(const double x,
									 const double y) {
		return 2.0 * maxwellTime_2(x,y) * (exp(constants.t/maxwellTime_2(x,y)) - 1.0) *
			(shearModulusRatio_1(x,y) * exp(constants.t * (maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_2(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_3(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y)))) *
			(constants.a * (2.0 * x - y) + constants.b * (2.0 * y - x)) *
			exp(-constants.t * (maxwellTime_1(x,y) * maxwellTime_2(x,y) + 2.0 * maxwellTime_1(x,y) * maxwellTime_3(x,y) + maxwellTime_2(x,y) * maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y) * maxwellTime_3(x,y)))/(3.0 * constants.t);
	} // viscousStrain_2_xx
	
	static double viscousStrain_2_yy(const double x,
									 const double y) {
		return -2.0 * maxwellTime_2(x,y) * (exp(constants.t/maxwellTime_2(x,y)) - 1.0) *
			(shearModulusRatio_1(x,y) * exp(constants.t * (maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_2(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_3(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y)))) *
			(constants.a * (x - 2.0 *y) - constants.b * (2.0 * x - y)) *
			exp(-constants.t * (maxwellTime_1(x,y) * maxwellTime_2(x,y) + 2.0 * maxwellTime_1(x,y) * maxwellTime_3(x,y) + maxwellTime_2(x,y) * maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y) * maxwellTime_3(x,y)))/(3.0 * constants.t);
	} // viscousStrain_2_yy
	
	static double viscousStrain_2_zz(const double x,
									 const double y) {
		return -2.0 * maxwellTime_2(x,y) * (exp(constants.t/maxwellTime_2(x,y)) - 1.0) *
			(shearModulusRatio_1(x,y) * exp(constants.t * (maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_2(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_3(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y)))) *
			(constants.a * (x + y) + constants.b * (x + y)) *
			exp(-constants.t * (maxwellTime_1(x,y) * maxwellTime_2(x,y) + 2.0 * maxwellTime_1(x,y) * maxwellTime_3(x,y) + maxwellTime_2(x,y) * maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y) * maxwellTime_3(x,y)))/(3.0 * constants.t);
	} // viscousStrain_2_zz
	
	static double viscousStrain_2_xy(const double x,
									 const double y) {
		return maxwellTime_2(x,y) * (exp(constants.t/maxwellTime_2(x,y)) - 1.0) *
			(shearModulusRatio_1(x,y) * exp(constants.t * (maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_2(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_3(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y)))) *
			(constants.b * (x + y) + constants.c * (x + y)) *
			exp(-constants.t * (maxwellTime_1(x,y) * maxwellTime_2(x,y) + 2.0 * maxwellTime_1(x,y) * maxwellTime_3(x,y) + maxwellTime_2(x,y) * maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y) * maxwellTime_3(x,y)))/constants.t;
	} // viscousStrain_2_xy
	
	// Viscous strain 3
	
	static double viscousStrain_3_xx(const double x,
									 const double y) {
		return 2.0 * maxwellTime_3(x,y) * (exp(constants.t/maxwellTime_3(x,y)) - 1.0) *
			(shearModulusRatio_1(x,y) * exp(constants.t * (maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_2(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_3(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y)))) *
			(constants.a * (2.0 * x - y) + constants.b * (2.0 * y - x)) *
			exp(-constants.t * (2.0 * maxwellTime_1(x,y) * maxwellTime_2(x,y) + maxwellTime_1(x,y) * maxwellTime_3(x,y) + maxwellTime_2(x,y) * maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y) * maxwellTime_3(x,y)))/(3.0 * constants.t);
	} // viscousStrain_3_xx
	
	static double viscousStrain_3_yy(const double x,
									 const double y) {
		return -2.0 * maxwellTime_3(x,y) * (exp(constants.t/maxwellTime_3(x,y)) - 1.0) *
			(shearModulusRatio_1(x,y) * exp(constants.t * (maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_2(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_3(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y)))) *
			(constants.a * (x - 2.0 *y) - constants.b * (2.0 * x - y)) *
			exp(-constants.t * (2.0 * maxwellTime_1(x,y) * maxwellTime_2(x,y) + maxwellTime_1(x,y) * maxwellTime_3(x,y) + maxwellTime_2(x,y) * maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y) * maxwellTime_3(x,y)))/(3.0 * constants.t);
	} // viscousStrain_3_yy
	
	static double viscousStrain_3_zz(const double x,
									 const double y) {
		return -2.0 * maxwellTime_3(x,y) * (exp(constants.t/maxwellTime_3(x,y)) - 1.0) *
			(shearModulusRatio_1(x,y) * exp(constants.t * (maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_2(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_3(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y)))) *
			(constants.a * (x + y) + constants.b * (x + y)) *
			exp(-constants.t * (2.0 * maxwellTime_1(x,y) * maxwellTime_2(x,y) + maxwellTime_1(x,y) * maxwellTime_3(x,y) + maxwellTime_2(x,y) * maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y) * maxwellTime_3(x,y)))/(3.0 * constants.t);
	} // viscousStrain_3_zz
	
	static double viscousStrain_3_xy(const double x,
									 const double y) {
		return maxwellTime_3(x,y) * (exp(constants.t/maxwellTime_3(x,y)) - 1.0) *
			(shearModulusRatio_1(x,y) * exp(constants.t * (maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_2(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_3(x,y))) +
			 shearModulusRatio_3(x,y) * exp(constants.t * (maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y)))) *
			(constants.b * (x + y) + constants.c * (x + y)) *
			exp(-constants.t * (2.0 * maxwellTime_1(x,y) * maxwellTime_2(x,y) + maxwellTime_1(x,y) * maxwellTime_3(x,y) + maxwellTime_2(x,y) * maxwellTime_3(x,y))/(maxwellTime_1(x,y) * maxwellTime_2(x,y) * maxwellTime_3(x,y)))/constants.t;
	} // viscousStrain_3_xy
	
	// Total strain for perturbed solution.
	
	static double totalStrainUpdate_xx(const double x,
									   const double y) {
		return (2.0*constants.a*x + 2.0*constants.b*y) *
			(shearModulusRatio_1(x,y) * exp(-constants.t/maxwellTime_1(x,y)) +
			 shearModulusRatio_2(x,y) * exp(-constants.t/maxwellTime_2(x,y)) +
			 shearModulusRatio_3(x,y) * exp(-constants.t/maxwellTime_3(x,y))) + constants.d;
	} // totalStrainUpdate_xx
	
	static double totalStrainUpdate_yy(const double x,
									   const double y) {
		return (2.0*constants.a*y + 2.0*constants.b*x) *
			(shearModulusRatio_1(x,y) * exp(-constants.t/maxwellTime_1(x,y)) +
			 shearModulusRatio_2(x,y) * exp(-constants.t/maxwellTime_2(x,y)) +
			 shearModulusRatio_3(x,y) * exp(-constants.t/maxwellTime_3(x,y)));
	} // totalStrainUpdate_yy
	
	static double totalStrainUpdate_zz(const double x,
									   const double y) {
		return 0.0;
	} // totalStrainUpdate_zz
	
	static double totalStrainUpdate_xy(const double x,
									   const double y) {
		return (shearModulusRatio_1(x,y)*exp(constants.t*(maxwellTime_2(x,y) + maxwellTime_3(x,y))/(maxwellTime_2(x,y)*maxwellTime_3(x,y))) +
				shearModulusRatio_2(x,y)*exp(constants.t*(maxwellTime_1(x,y) + maxwellTime_3(x,y))/(maxwellTime_1(x,y)*maxwellTime_3(x,y))) +
				shearModulusRatio_3(x,y)*exp(constants.t*(maxwellTime_1(x,y) + maxwellTime_2(x,y))/(maxwellTime_1(x,y)*maxwellTime_2(x,y))))*
			(constants.b * (x + y) + constants.c * (x + y)) *
			exp(-constants.t*(maxwellTime_1(x,y)*maxwellTime_2(x,y) + maxwellTime_1(x,y)*maxwellTime_3(x,y) + maxwellTime_2(x,y)*maxwellTime_3(x,y))/(maxwellTime_1(x,y)*maxwellTime_2(x,y)*maxwellTime_3(x,y))) + constants.d/2.0;
	} // totalStrainUpdate_xy

	// Values needed to compute viscous strain for perturbed solution.
	static double meanStrainT(const double x,
							  const double y) {
		return (totalStrain_xx(x,y) + totalStrain_yy(x,y))/3.0;
	}

	static double devStrainT_xx(const double x,
								const double y) {
		return totalStrain_xx(x,y) - meanStrainT(x,y);
	}

	static double devStrainT_yy(const double x,
								const double y) {
		return totalStrain_yy(x,y) - meanStrainT(x,y);
	}

	static double devStrainT_zz(const double x,
								const double y) {
		return totalStrain_zz(x,y) - meanStrainT(x,y);
	}

	static double devStrainT_xy(const double x,
								const double y) {
		return totalStrain_xy(x,y);
	}

	static double meanStrainTplusDt(const double x,
									const double y) {
		return (totalStrainUpdate_xx(x,y) + totalStrainUpdate_yy(x,y))/3.0;
	}

	static double devStrainTplusDt_xx(const double x,
									  const double y) {
		return totalStrainUpdate_xx(x,y) - meanStrainTplusDt(x,y);
	}

	static double devStrainTplusDt_yy(const double x,
									  const double y) {
		return totalStrainUpdate_yy(x,y) - meanStrainTplusDt(x,y);
	}

	static double devStrainTplusDt_zz(const double x,
									  const double y) {
		return totalStrainUpdate_zz(x,y) - meanStrainTplusDt(x,y);
	}

	static double devStrainTplusDt_xy(const double x,
								const double y) {
		return totalStrain_xy(x,y);
	}

	static double dq_1(const double x,
					   const double y) {
		return maxwellTime_1(x,y) * (1.0 - exp(-constants.dt/maxwellTime_1(x,y)))/constants.dt;
	}

	static double dq_2(const double x,
					   const double y) {
		return maxwellTime_2(x,y) * (1.0 - exp(-constants.dt/maxwellTime_2(x,y)))/constants.dt;
	}

	static double dq_3(const double x,
					   const double y) {
		return maxwellTime_3(x,y) * (1.0 - exp(-constants.dt/maxwellTime_3(x,y)))/constants.dt;
	}

	// Viscous strain for perturbed solution.

    static double viscousStrainUpdate_1_xx(const double x,
										   const double y) {
		return exp(-constants.dt/maxwellTime_1(x,y)) * viscousStrain_1_xx(x,y) + dq_1(x,y) * (devStrainTplusDt_xx(x,y) - devStrainT_xx(x,y));
    } // viscousStrainUpdate_1_xx

    static double viscousStrainUpdate_1_yy(const double x,
										   const double y) {
		return exp(-constants.dt/maxwellTime_1(x,y)) * viscousStrain_1_yy(x,y) + dq_1(x,y) * (devStrainTplusDt_yy(x,y) - devStrainT_yy(x,y));
    } // viscousStrainUpdate_1_yy

    static double viscousStrainUpdate_1_zz(const double x,
										   const double y) {
		return exp(-constants.dt/maxwellTime_1(x,y)) * viscousStrain_1_zz(x,y) + dq_1(x,y) * (devStrainTplusDt_zz(x,y) - devStrainT_zz(x,y));
    } // viscousStrainUpdate_1_zz

    static double viscousStrainUpdate_1_xy(const double x,
										   const double y) {
		return exp(-constants.dt/maxwellTime_1(x,y)) * viscousStrain_1_xy(x,y) + dq_1(x,y) * (devStrainTplusDt_xy(x,y) - devStrainT_xy(x,y));
    } // viscousStrainUpdate_1_xy

    static double viscousStrainUpdate_2_xx(const double x,
										   const double y) {
		return exp(-constants.dt/maxwellTime_2(x,y)) * viscousStrain_2_xx(x,y) + dq_2(x,y) * (devStrainTplusDt_xx(x,y) - devStrainT_xx(x,y));
    } // viscousStrainUpdate_2_xx

    static double viscousStrainUpdate_2_yy(const double x,
										   const double y) {
		return exp(-constants.dt/maxwellTime_2(x,y)) * viscousStrain_2_yy(x,y) + dq_2(x,y) * (devStrainTplusDt_yy(x,y) - devStrainT_yy(x,y));
    } // viscousStrainUpdate_2_yy

    static double viscousStrainUpdate_2_zz(const double x,
										   const double y) {
		return exp(-constants.dt/maxwellTime_2(x,y)) * viscousStrain_2_zz(x,y) + dq_2(x,y) * (devStrainTplusDt_zz(x,y) - devStrainT_zz(x,y));
    } // viscousStrainUpdate_2_zz

    static double viscousStrainUpdate_2_xy(const double x,
										   const double y) {
		return exp(-constants.dt/maxwellTime_2(x,y)) * viscousStrain_2_xy(x,y) + dq_2(x,y) * (devStrainTplusDt_xy(x,y) - devStrainT_xy(x,y));
    } // viscousStrainUpdate_2_xy

    static double viscousStrainUpdate_3_xx(const double x,
										   const double y) {
		return exp(-constants.dt/maxwellTime_3(x,y)) * viscousStrain_3_xx(x,y) + dq_3(x,y) * (devStrainTplusDt_xx(x,y) - devStrainT_xx(x,y));
    } // viscousStrainUpdate_3_xx

    static double viscousStrainUpdate_3_yy(const double x,
										   const double y) {
		return exp(-constants.dt/maxwellTime_3(x,y)) * viscousStrain_3_yy(x,y) + dq_3(x,y) * (devStrainTplusDt_yy(x,y) - devStrainT_yy(x,y));
    } // viscousStrainUpdate_3_yy

    static double viscousStrainUpdate_3_zz(const double x,
										   const double y) {
		return exp(-constants.dt/maxwellTime_3(x,y)) * viscousStrain_3_zz(x,y) + dq_3(x,y) * (devStrainTplusDt_zz(x,y) - devStrainT_zz(x,y));
    } // viscousStrainUpdate_3_zz

    static double viscousStrainUpdate_3_xy(const double x,
										   const double y) {
		return exp(-constants.dt/maxwellTime_3(x,y)) * viscousStrain_3_xy(x,y) + dq_3(x,y) * (devStrainTplusDt_xy(x,y) - devStrainT_xy(x,y));
    } // viscousStrainUpdate_3_xy
	
	// Body force
	static double bodyforce_x(const double x,
							  const double y) {
		return 6.0*(constants.a + constants.b)*bulkModulus(x,y)*shearModulusRatio_1(x,y)*exp(-constants.t/maxwellTime_1(x,y)) +
			6.0*(constants.a + constants.b)*bulkModulus(x,y)*shearModulusRatio_2(x,y)*exp(-constants.t/maxwellTime_2(x,y)) +
			6.0*(constants.a + constants.b)*bulkModulus(x,y)*shearModulusRatio_3(x,y)*exp(-constants.t/maxwellTime_3(x,y));
	} // bodyforce_x
	
	static double bodyforce_y(const double x,
							  const double y) {
		return 6.0*(constants.a + constants.b)*bulkModulus(x,y)*shearModulusRatio_1(x,y)*exp(-constants.t/maxwellTime_1(x,y)) +
			6.0*(constants.a + constants.b)*bulkModulus(x,y)*shearModulusRatio_2(x,y)*exp(-constants.t/maxwellTime_2(x,y)) +
			6.0*(constants.a + constants.b)*bulkModulus(x,y)*shearModulusRatio_3(x,y)*exp(-constants.t/maxwellTime_3(x,y));
	} // bodyforce_y
	
	static const char* bodyforce_units(void) {
		return "kg*m/s**2";
	}
	
	// Spatial database user functions for solution subfields.
	
	// Displacement
	static double disp_x(const double x,
						 const double y) {
		return (constants.a*x*x + 2.0*constants.b*x*y + constants.c*y*y) *
			(shearModulusRatio_1(x,y)*exp(-constants.t/maxwellTime_1(x,y)) +
			 shearModulusRatio_2(x,y)*exp(-constants.t/maxwellTime_2(x,y)) +
			 shearModulusRatio_3(x,y)*exp(-constants.t/maxwellTime_3(x,y)));
	} // disp_x
	
	static double disp_y(const double x,
						 const double y) {
		return (constants.a*y*y + 2.0*constants.b*x*y + constants.c*x*x) *
			(shearModulusRatio_1(x,y)*exp(-constants.t/maxwellTime_1(x,y)) +
			 shearModulusRatio_2(x,y)*exp(-constants.t/maxwellTime_2(x,y)) +
			 shearModulusRatio_3(x,y)*exp(-constants.t/maxwellTime_3(x,y)));
	} // disp_y
	
	static const char* disp_units(void) {
		return "m";
	} // disp_units
	
	static double disp_dot_x(const double x,
							 const double y) {
		return -(constants.a*x*x + 2.0*constants.b*x*y + constants.c*y*y) *
			(shearModulusRatio_1(x,y)*exp(-constants.t/maxwellTime_1(x,y))/maxwellTime_1(x,y) +
			 shearModulusRatio_2(x,y)*exp(-constants.t/maxwellTime_2(x,y))/maxwellTime_2(x,y) +
			 shearModulusRatio_3(x,y)*exp(-constants.t/maxwellTime_3(x,y))/maxwellTime_3(x,y));
	} // disp_dot_x
	static double disp_dot_y(const double x,
							 const double y) {
		return -(constants.a*y*y + 2.0*constants.b*x*y + constants.c*x*x) *
			(shearModulusRatio_1(x,y)*exp(-constants.t/maxwellTime_1(x,y))/maxwellTime_1(x,y) +
			 shearModulusRatio_2(x,y)*exp(-constants.t/maxwellTime_2(x,y))/maxwellTime_2(x,y) +
			 shearModulusRatio_3(x,y)*exp(-constants.t/maxwellTime_3(x,y))/maxwellTime_3(x,y));
	} // disp_dot_y
	
	static const char* disp_dot_units(void) {
		return "m/s";
	} // disp_dot_units
	
	// Displacement + perturbation
	static double disp_perturb_x(const double x,
								 const double y) {
		return disp_x(x, y) + constants.d * x;
	} // disp_perturb_x
	
	static double disp_perturb_y(const double x,
								 const double y) {
		return disp_y(x, y) + constants.d * x;
	} // disp_perturb_y
	
protected:
	void setUp(void) {
		TestIsotropicLinearGenMaxwellPlaneStrain::setUp();
		_mydata = new TestIsotropicLinearGenMaxwellPlaneStrain_Data(); CPPUNIT_ASSERT(_mydata);
		
		// dimension set in base class.
		// meshFilename set in derived class.
		_mydata->boundaryLabel = "boundary";
		
		CPPUNIT_ASSERT(_mydata->normalizer);
		_mydata->normalizer->setLengthScale(1.0e+03);
		_mydata->normalizer->setTimeScale(2.0e+7);
		_mydata->normalizer->setDensityScale(3.0e+3);
		_mydata->normalizer->setPressureScale(2.25e+10);
		
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
		_mydata->auxDB->addValue("viscous_strain_2_xx", viscousStrain_2_xx, "none");
		_mydata->auxDB->addValue("viscous_strain_2_yy", viscousStrain_2_yy, "none");
		_mydata->auxDB->addValue("viscous_strain_2_zz", viscousStrain_2_zz, "none");
		_mydata->auxDB->addValue("viscous_strain_2_xy", viscousStrain_2_xy, "none");
		_mydata->auxDB->addValue("viscous_strain_3_xx", viscousStrain_3_xx, "none");
		_mydata->auxDB->addValue("viscous_strain_3_yy", viscousStrain_3_yy, "none");
		_mydata->auxDB->addValue("viscous_strain_3_zz", viscousStrain_3_zz, "none");
		_mydata->auxDB->addValue("viscous_strain_3_xy", viscousStrain_3_xy, "none");
		_mydata->auxDB->addValue("total_strain_xx", totalStrain_xx, "none");
		_mydata->auxDB->addValue("total_strain_yy", totalStrain_yy, "none");
		_mydata->auxDB->addValue("total_strain_zz", totalStrain_zz, "none");
		_mydata->auxDB->addValue("total_strain_xy", totalStrain_xy, "none");
		_mydata->auxDB->addValue("body_force_x", bodyforce_x, bodyforce_units());
		_mydata->auxDB->addValue("body_force_y", bodyforce_y, bodyforce_units());

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
		_mydata->auxUpdateDB->addValue("viscous_strain_2_xx", viscousStrainUpdate_2_xx, "none");
		_mydata->auxUpdateDB->addValue("viscous_strain_2_yy", viscousStrainUpdate_2_yy, "none");
		_mydata->auxUpdateDB->addValue("viscous_strain_2_zz", viscousStrainUpdate_2_zz, "none");
		_mydata->auxUpdateDB->addValue("viscous_strain_2_xy", viscousStrainUpdate_2_xy, "none");
		_mydata->auxUpdateDB->addValue("viscous_strain_3_xx", viscousStrainUpdate_3_xx, "none");
		_mydata->auxUpdateDB->addValue("viscous_strain_3_yy", viscousStrainUpdate_3_yy, "none");
		_mydata->auxUpdateDB->addValue("viscous_strain_3_zz", viscousStrainUpdate_3_zz, "none");
		_mydata->auxUpdateDB->addValue("viscous_strain_3_xy", viscousStrainUpdate_3_xy, "none");
		_mydata->auxUpdateDB->addValue("total_strain_xx", totalStrainUpdate_xx, "none");
		_mydata->auxUpdateDB->addValue("total_strain_yy", totalStrainUpdate_yy, "none");
		_mydata->auxUpdateDB->addValue("total_strain_zz", totalStrainUpdate_zz, "none");
		_mydata->auxUpdateDB->addValue("total_strain_xy", totalStrainUpdate_xy, "none");
		_mydata->auxUpdateDB->addValue("body_force_x", bodyforce_x, bodyforce_units());
		_mydata->auxUpdateDB->addValue("body_force_y", bodyforce_y, bodyforce_units());
		
		CPPUNIT_ASSERT(_mydata->solnDB);
		_mydata->solnDB->addValue("displacement_x", disp_x, disp_units());
		_mydata->solnDB->addValue("displacement_y", disp_y, disp_units());
		_mydata->solnDB->addValue("displacement_dot_x", disp_dot_x, disp_dot_units());
		_mydata->solnDB->addValue("displacement_dot_y", disp_dot_y, disp_dot_units());
		
		CPPUNIT_ASSERT(_mydata->perturbDB);
		_mydata->perturbDB->addValue("displacement_x", disp_perturb_x, disp_units());
		_mydata->perturbDB->addValue("displacement_y", disp_perturb_y, disp_units());
		_mydata->perturbDB->addValue("displacement_dot_x", disp_dot_x, disp_dot_units());
		_mydata->perturbDB->addValue("displacement_dot_y", disp_dot_y, disp_dot_units());
		
		CPPUNIT_ASSERT(_mymaterial);
		_mymaterial->useInertia(false);
		_mymaterial->useBodyForce(true);
		_mymaterial->useReferenceState(false);
		
		_mymaterial->setLabel("Isotropic Linear Generalized Maxwell Plane Strain");
		_mymaterial->id(24);
		
	} // setUp

}; // TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain
const double pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain::SMALL = 1.0e-5;

const pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain::AuxConstants pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain::constants = {
	1.0e-7, // a
	2.5e-7, // b
	3.0e-7, // c
	9.0e-8, // d
	5.0e+7, // t
	5.0e+7, // dt
};


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP1 :
	public pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain {

	CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP1,
						   TestIsotropicLinearGenMaxwellPlaneStrain);
	CPPUNIT_TEST_SUITE_END();

	void setUp(void) {
		TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain::setUp();
		CPPUNIT_ASSERT(_mydata);

		_mydata->meshFilename = "data/tri_small.mesh";

		_mydata->numSolnSubfields = 1;
		static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
			pylith::topology::Field::Discretization(1, 1), // disp
		};
		_mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

		_initializeMin();
	} // setUp

}; // TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP1
// Do not use this test because TriP1 cannot represent a linear strain field.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP1);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP2 :
	public pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain {

	CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP2,
						   TestIsotropicLinearGenMaxwellPlaneStrain);
	CPPUNIT_TEST_SUITE_END();

	void setUp(void) {
		TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain::setUp();
		CPPUNIT_ASSERT(_mydata);

		_mydata->meshFilename = "data/tri_small.mesh";

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

}; // TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP2);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP3 :
	public pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain {

	CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP3,
						   TestIsotropicLinearGenMaxwellPlaneStrain);
	CPPUNIT_TEST_SUITE_END();

	void setUp(void) {
		TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain::setUp();
		CPPUNIT_ASSERT(_mydata);

		_mydata->meshFilename = "data/tri_small.mesh";

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

}; // TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP3
// Remove this test for now until higher order integration is done properly.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP3);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP4 :
	public pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain {

	CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP4,
						   TestIsotropicLinearGenMaxwellPlaneStrain);
	CPPUNIT_TEST_SUITE_END();

	void setUp(void) {
		TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain::setUp();
		CPPUNIT_ASSERT(_mydata);

		_mydata->meshFilename = "data/tri_small.mesh";

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

}; // TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP4
// Remove this test for now until higher order integration is done properly.
//CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_TriP4);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ1 :
	public pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain {

	CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ1,
						   TestIsotropicLinearGenMaxwellPlaneStrain);
	CPPUNIT_TEST_SUITE_END();

	void setUp(void) {
		TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain::setUp();
		CPPUNIT_ASSERT(_mydata);

		_mydata->meshFilename = "data/quad_aligned.mesh";

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

}; // TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ1
// Do not use this test because TriP1 cannot represent a linear strain field.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ1);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ2 :
	public pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain {

	CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ2,
						   TestIsotropicLinearGenMaxwellPlaneStrain);
	CPPUNIT_TEST_SUITE_END();

	void setUp(void) {
		TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain::setUp();
		CPPUNIT_ASSERT(_mydata);

		_mydata->meshFilename = "data/quad_aligned.mesh";

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

}; // TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ2);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ3 :
	public pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain {

	CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ3,
						   TestIsotropicLinearGenMaxwellPlaneStrain);
	CPPUNIT_TEST_SUITE_END();

	void setUp(void) {
		TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain::setUp();
		CPPUNIT_ASSERT(_mydata);

		_mydata->meshFilename = "data/quad_aligned.mesh";

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

}; // TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ3
// Remove this test for now until higher order integration is done properly.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ3);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ4 :
	public pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain {

	CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ4,
						   TestIsotropicLinearGenMaxwellPlaneStrain);
	CPPUNIT_TEST_SUITE_END();

	void setUp(void) {
		TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain::setUp();
		CPPUNIT_ASSERT(_mydata);

		_mydata->meshFilename = "data/quad_aligned.mesh";

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

}; // TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ4
// Remove this test for now until higher order integration is done properly.
//CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain_QuadQ4);


// End of file
