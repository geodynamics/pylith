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

#include "TestIsotropicLinearMaxwellPlaneStrain.hh" // Implementation of cases

#include "pylith/materials/IsotropicLinearMaxwellPlaneStrain.hh" // USES IsotropicLinearMaxwellPlaneStrain
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // :TEMPORARY: USES PYLITH_JOURNAL_ERROR

// forward declarations
namespace pylith {
    namespace materials {
        class TestIsotropicLinearMaxwellPlaneStrain_LinearStrain;

        class TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP1;
        class TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP2;
        class TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP3;
        class TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP4;

        class TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ1;
        class TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ2;
        class TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ3;
        class TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ4;

    } // materials
} // pylith


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain :
    public pylith::materials::TestIsotropicLinearMaxwellPlaneStrain {

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

    // Viscosity
    static double viscosity(const double x,
                            const double y) {
        return 7.91700159488e+19;
    } // viscosity
    static const char* viscosity_units(void) {
        return "Pa*s";
    } // viscosity_units

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

    // Maxwell time
    static double maxwellTime(const double x,
                              const double y) {
        return viscosity(x,y) / shearModulus(x,y);
    } // maxwellTime
    static const char* maxwellTime_units(void) {
        return "s";
    } // maxwellTime_units

    // Temporal and spatial constants.
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
        return (2.0*constants.a*x + 2.0*constants.b*y) * exp(-constants.t/maxwellTime(x,y));
    } // totalStrain_xx

    static double totalStrain_yy(const double x,
                                 const double y) {
        return (2.0*constants.b*x + 2.0*constants.a*y) * exp(-constants.t/maxwellTime(x,y));
    } // totalStrain_yy

    static double totalStrain_zz(const double x,
                                 const double y) {
        return 0.0;
    } // totalStrain_zz

    static double totalStrain_xy(const double x,
                                 const double y) {
        return (constants.b*(x+y) + constants.c*(x+y))*exp(-constants.t/maxwellTime(x,y));
    } // totalStrain_xy

    // Viscous strain

    static double viscousStrain_xx(const double x,
                                   const double y) {
        return 2.0*maxwellTime(x,y)*(exp(constants.t/maxwellTime(x,y)) - 1.0)*(constants.a*(2.0*x-y) - constants.b*(x-2.0*y)) * exp(-2.0*constants.t/maxwellTime(x,y))/(3.0*constants.t);
    } // viscousStrain_xx

    static double viscousStrain_yy(const double x,
                                   const double y) {
        return -2.0*maxwellTime(x,y)*(exp(constants.t/maxwellTime(x,y)) - 1.0)*(constants.a*(x-2.0*y) - constants.b*(2.0*x-y)) * exp(-2.0*constants.t/maxwellTime(x,y))/(3.0*constants.t);
    } // viscousStrain_yy

    static double viscousStrain_zz(const double x,
                                   const double y) {
        return -2.0*maxwellTime(x,y)*(exp(constants.t/maxwellTime(x,y)) - 1.0)*(constants.a*(x+y) + constants.b*(x+y)) * exp(-2.0*constants.t/maxwellTime(x,y))/(3.0*constants.t);
    } // viscousStrain_zz

    static double viscousStrain_xy(const double x,
                                   const double y) {
        return maxwellTime(x,y)*(exp(constants.t/maxwellTime(x,y)) - 1.0)*(constants.b*(x+y) + constants.c*(x+y)) * exp(-2.0*constants.t/maxwellTime(x,y))/constants.t;
    } // viscousStrain_xy

    // Total strain for perturbed solution.

    static double totalStrainUpdate_xx(const double x,
									   const double y) {
        return (2.0*constants.a*x + 2.0*constants.b*y) * exp(-constants.t/maxwellTime(x,y)) + constants.d;
    } // totalStrainUpdate_xx

    static double totalStrainUpdate_yy(const double x,
									   const double y) {
        return (2.0*constants.b*x + 2.0*constants.a*y) * exp(-constants.t/maxwellTime(x,y));
    } // totalStrainUpdate_yy

    static double totalStrainUpdate_zz(const double x,
									   const double y) {
        return 0.0;
    } // totalStrainUpdate_zz

    static double totalStrainUpdate_xy(const double x,
									   const double y) {
        return (constants.b*(x+y) + constants.c*(x+y))*exp(-constants.t/maxwellTime(x,y)) + constants.d/2.0;
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
		return totalStrainUpdate_xy(x,y);
	}

	static double dq(const double x,
					 const double y) {
		return maxwellTime(x,y) * (1.0 - exp(-constants.dt/maxwellTime(x,y)))/constants.dt;
	}
	
    // Viscous strain for perturbed solution.

    static double viscousStrainUpdate_xx(const double x,
										 const double y) {
		return exp(-constants.dt/maxwellTime(x,y)) * viscousStrain_xx(x,y) + dq(x,y) * (devStrainTplusDt_xx(x,y) - devStrainT_xx(x,y));
    } // viscousStrainUpdate_xx

    static double viscousStrainUpdate_yy(const double x,
										 const double y) {
		return exp(-constants.dt/maxwellTime(x,y)) * viscousStrain_yy(x,y) + dq(x,y) * (devStrainTplusDt_yy(x,y) - devStrainT_yy(x,y));
    } // viscousStrainUpdate_yy

    static double viscousStrainUpdate_zz(const double x,
										 const double y) {
		return exp(-constants.dt/maxwellTime(x,y)) * viscousStrain_zz(x,y) + dq(x,y) * (devStrainTplusDt_zz(x,y) - devStrainT_zz(x,y));
    } // viscousStrainUpdate_zz

    static double viscousStrainUpdate_xy(const double x,
										 const double y) {
		return exp(-constants.dt/maxwellTime(x,y)) * viscousStrain_xy(x,y) + dq(x,y) * (devStrainTplusDt_xy(x,y) - devStrainT_xy(x,y));
    } // viscousStrainUpdate_xy

/*
    static double viscousStrainUpdate_xx(const double x,
										 const double y) {
		return 2.0*maxwellTime(x,y)*(constants.d*exp((constants.dt + 2.0*constants.t)/maxwellTime(x,y)) + (2.0*constants.a*x - constants.a*y - constants.b*x + 2.0*constants.b*y)*exp(constants.t/maxwellTime(x,y)))*(exp(constants.t/maxwellTime(x,y)) - 1.0)*exp((-constants.dt - 3.0*constants.t)/maxwellTime(x,y))/(3.0*constants.t);
    } // viscousStrainUpdate_xx

    static double viscousStrainUpdate_yy(const double x,
										 const double y) {
		return -maxwellTime(x,y)*(constants.d*exp((constants.dt + 2.0*constants.t)/maxwellTime(x,y)) + 2.0*(constants.a*x - 2.0*constants.a*y - 2.0*constants.b*x + constants.b*y)*exp(constants.t/maxwellTime(x,y)))*(exp(constants.t/maxwellTime(x,y)) - 1.0)*exp(-(constants.dt + 3.0*constants.t)/maxwellTime(x,y))/(3.0*constants.t);
    } // viscousStrainUpdate_yy

    static double viscousStrainUpdate_zz(const double x,
										 const double y) {
		return -maxwellTime(x,y)*(constants.d*exp((constants.dt + 2.0*constants.t)/maxwellTime(x,y)) + 2.0*(constants.a*x + constants.a*y + constants.b*x + constants.b*y)*exp(constants.t/maxwellTime(x,y)))*(exp(constants.t/maxwellTime(x,y)) - 1.0)*exp(-(constants.dt + 3.0*constants.t)/maxwellTime(x,y))/(3.0*constants.t);
    } // viscousStrainUpdate_zz

    static double viscousStrainUpdate_xy(const double x,
										 const double y) {
		return maxwellTime(x,y)*(constants.d*exp((constants.dt + 2.0*constants.t)/maxwellTime(x,y)) + 2.0*(constants.b*x + constants.b*y + constants.c*x + constants.c*y)*exp(constants.t/maxwellTime(x,y)))*(exp(constants.t/maxwellTime(x,y)) - 1.0)*exp((-constants.dt - 3.0*constants.t)/maxwellTime(x,y))/(2.0*constants.t);
    } // viscousStrainUpdate_xy
*/
    // Body force
    static double bodyforce_x(const double x,
                              const double y) {
        return 6.0*bulkModulus(x,y)*(constants.a + constants.b) * exp(-constants.t/maxwellTime(x,y));
    } // bodyforce_x

    static double bodyforce_y(const double x,
                              const double y) {
        return 6.0*bulkModulus(x,y)*(constants.a + constants.b) * exp(-constants.t/maxwellTime(x,y));
    } // bodyforce_y

    static const char* bodyforce_units(void) {
        return "kg*m/s**2";
    }

    // Spatial database user functions for solution subfields.

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return (constants.a*x*x + 2.0*constants.b*x*y + constants.c*y*y) * exp(-constants.t/maxwellTime(x,y));
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return (constants.a*y*y + 2.0*constants.b*x*y + constants.c*x*x) * exp(-constants.t/maxwellTime(x,y));
    } // disp_y

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    static double disp_dot_x(const double x,
                             const double y) {
        return -(constants.a*x*x + 2.0*constants.b*x*y + constants.c*y*y) * exp(-constants.t/maxwellTime(x,y)) / maxwellTime(x,y);
    } // disp_dot_x

    static double disp_dot_y(const double x,
                             const double y) {
        return -(constants.a*y*y + 2.0*constants.b*x*y + constants.c*x*x) * exp(-constants.t/maxwellTime(x,y))/maxwellTime(x,y);
    } // disp_dot_y

    static const char* disp_dot_units(void) {
        return "m/s";
    } // disp_dot_units

    // Displacement + perturbation
    static double disp_perturb_x(const double x,
                                 const double y) {
        return disp_x(x, y) + constants.d*x;
    } // disp_perturb_x

    static double disp_perturb_y(const double x,
                                 const double y) {
        return disp_y(x, y) + constants.d*x;
    } // disp_perturb_y

protected:
    void setUp(void) {
        TestIsotropicLinearMaxwellPlaneStrain::setUp();
        _mydata = new TestIsotropicLinearMaxwellPlaneStrain_Data(); CPPUNIT_ASSERT(_mydata);

        // dimension set in base class.
        // meshFilename set in derived class.
        _mydata->boundaryLabel = "boundary";

        CPPUNIT_ASSERT(_mydata->normalizer);
        _mydata->normalizer->setLengthScale(1.0e+03);
        _mydata->normalizer->setTimeScale(6.3e+8);
        _mydata->normalizer->setDensityScale(4.0e+3);
        _mydata->normalizer->setPressureScale(2.5e+11);

        _mydata->t = constants.t/_mydata->normalizer->getTimeScale();
        _mydata->dt = constants.dt/_mydata->normalizer->getTimeScale();
        _mydata->s_tshift = 1.0 / _mydata->dt;

        // solnDiscretizations set in derived class.

        _mydata->numAuxSubfields = 7;
        static const char* _auxSubfields[7] =
			{"density", "shear_modulus", "bulk_modulus", "maxwell_time",
			 "viscous_strain", "total_strain", "body_force"};
        _mydata->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // maxwell_time
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
        _mydata->auxDB->addValue("viscosity", viscosity, viscosity_units());
        _mydata->auxDB->addValue("maxwell_time", maxwellTime, maxwellTime_units());
        _mydata->auxDB->addValue("viscous_strain_xx", viscousStrain_xx, "none");
        _mydata->auxDB->addValue("viscous_strain_yy", viscousStrain_yy, "none");
        _mydata->auxDB->addValue("viscous_strain_zz", viscousStrain_zz, "none");
        _mydata->auxDB->addValue("viscous_strain_xy", viscousStrain_xy, "none");
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
        _mydata->auxUpdateDB->addValue("viscosity", viscosity, viscosity_units());
        _mydata->auxUpdateDB->addValue("maxwell_time", maxwellTime, maxwellTime_units());
        _mydata->auxUpdateDB->addValue("viscous_strain_xx", viscousStrainUpdate_xx, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_yy", viscousStrainUpdate_yy, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_zz", viscousStrainUpdate_zz, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_xy", viscousStrainUpdate_xy, "none");
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

        _mymaterial->setLabel("Isotropic Linear Maxwell Plane Strain");
        _mymaterial->id(24);

    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain
const double pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::SMALL = 1.0e-5;

const pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::AuxConstants pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::constants = {
    1.0e-7, // a
    2.5e-7, // b
    3.0e-7, // c
    9.0e-8, // d
    9.0e+7, // t
    5.0e+7  // dt
};


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP1 :
    public pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP1,
                           TestIsotropicLinearMaxwellPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tri_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP1
// Do not use this test because TriP1 cannot represent a linear strain field.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP1);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP2 :
    public pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP2,
                           TestIsotropicLinearMaxwellPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tri_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // maxwell_time
            pylith::topology::Field::Discretization(1, 2, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 2, false), // total_strain
            pylith::topology::Field::Discretization(0, 2), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP2);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP3 :
    public pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP3,
                           pylith::materials::TestIsotropicLinearMaxwellPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tri_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // maxwell_time
            pylith::topology::Field::Discretization(1, 3, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 3, false), // total_strain
            pylith::topology::Field::Discretization(0, 3), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP3
// Remove this test for now until higher order integration is done properly.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP3);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP4 :
    public pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP4,
                           TestIsotropicLinearMaxwellPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tri_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // maxwell_time
            pylith::topology::Field::Discretization(1, 4, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 4, false), // total_strain
            pylith::topology::Field::Discretization(0, 4), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP4
// Remove this test for now until higher order integration is done properly.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP4);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ1 :
    public pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ1,
                           TestIsotropicLinearMaxwellPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/quad_aligned.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // maxwell_time
            pylith::topology::Field::Discretization(1, 1, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 1, false), // total_strain
            pylith::topology::Field::Discretization(0, 1), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ1
// Do not use this test because TriP1 cannot represent a linear strain field.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ1);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ2 :
    public pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ2,
                           TestIsotropicLinearMaxwellPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/quad_aligned.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // maxwell_time
            pylith::topology::Field::Discretization(1, 2, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 2, false), // total_strain
            pylith::topology::Field::Discretization(0, 2), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ2);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ3 :
    public pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ3,
                           TestIsotropicLinearMaxwellPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/quad_aligned.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // maxwell_time
            pylith::topology::Field::Discretization(1, 3, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 3, false), // total_strain
            pylith::topology::Field::Discretization(0, 3), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ3
// Remove this test for now until higher order integration is done properly.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ3);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ4 :
    public pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ4,
                           TestIsotropicLinearMaxwellPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/quad_aligned.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // maxwell_time
            pylith::topology::Field::Discretization(1, 4, false), // viscous_strain
            pylith::topology::Field::Discretization(1, 4, false), // total_strain
            pylith::topology::Field::Discretization(0, 4), // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ4
// Remove this test for now until higher order integration is done properly.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ4);


// End of file
