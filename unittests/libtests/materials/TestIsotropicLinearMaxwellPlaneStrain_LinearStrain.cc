// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
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
        return 2.0*(exp(constants.t/maxwellTime(x,y)) - 1.0)*(constants.a*(2.0*x-y) - constants.b*(x-2.0*y)) * exp(-2.0*constants.t/maxwellTime(x,y))/3.0;
    } // viscousStrain_xx
    static double viscousStrain_yy(const double x,
                                   const double y) {
        return -2.0*(exp(constants.t/maxwellTime(x,y)) - 1.0)*(constants.a*(x-2.0*y) - constants.b*(2.0*x-y)) * exp(-2.0*constants.t/maxwellTime(x,y))/3.0;
    } // viscousStrain_yy
    static double viscousStrain_zz(const double x,
                                   const double y) {
        return -2.0*(exp(constants.t/maxwellTime(x,y)) - 1.0)*(constants.a*(x+y) + constants.b*(x+y)) * exp(-2.0*constants.t/maxwellTime(x,y))/3.0;
    } // viscousStrain_zz
    static double viscousStrain_xy(const double x,
                                   const double y) {
        return (exp(constants.t/maxwellTime(x,y)) - 1.0)*(constants.b*(x+y) + constants.c*(x+y)) * exp(-2.0*constants.t/maxwellTime(x,y));
    } // viscousStrain_xy

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
        return disp_x(x, y) + SMALL;
    } // disp_perturb_x
    static double disp_perturb_y(const double x,
                                 const double y) {
        return disp_y(x, y) + SMALL;
    } // disp_perturb_

protected:
    void setUp(void) {
        TestIsotropicLinearMaxwellPlaneStrain::setUp();
        _mydata = new TestIsotropicLinearMaxwellPlaneStrain_Data(); CPPUNIT_ASSERT(_mydata);

        // dimension set in base class.
        // meshFilename set in derived class.
        _mydata->boundaryLabel = "boundary";

        CPPUNIT_ASSERT(_mydata->normalizer);
        _mydata->normalizer->lengthScale(1.0e+03);
        _mydata->normalizer->timeScale(6.3e+8);
        _mydata->normalizer->densityScale(4.0e+3);
        _mydata->normalizer->pressureScale(2.5e+11);

        _mydata->t = constants.t;
        _mydata->dt = constants.dt;
        _mydata->tshift = 1.0 / _mydata->dt;

        // solnDiscretizations set in derived class.

        _mydata->numAuxSubfields = 7;
        static const char* _auxSubfields[7] = {"density", "shear_modulus", "bulk_modulus", "maxwell_time",
                                               "total_strain", "viscous_strain", "body_force"};
        _mydata->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
            {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
            {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
            {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // maxwell_time
            {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // total_strain
            {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // viscous_strain
            {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // body_force
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
        _mydata->auxDB->addValue("total_strain_xx", totalStrain_xx, "none");
        _mydata->auxDB->addValue("total_strain_yy", totalStrain_yy, "none");
        _mydata->auxDB->addValue("total_strain_zz", totalStrain_zz, "none");
        _mydata->auxDB->addValue("total_strain_xy", totalStrain_xy, "none");
        _mydata->auxDB->addValue("viscous_strain_xx", viscousStrain_xx, "none");
        _mydata->auxDB->addValue("viscous_strain_yy", viscousStrain_yy, "none");
        _mydata->auxDB->addValue("viscous_strain_zz", viscousStrain_zz, "none");
        _mydata->auxDB->addValue("viscous_strain_xy", viscousStrain_xy, "none");
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
        _mydata->auxUpdateDB->addValue("total_strain_xx", totalStrain_xx, "none");
        _mydata->auxUpdateDB->addValue("total_strain_yy", totalStrain_yy, "none");
        _mydata->auxUpdateDB->addValue("total_strain_zz", totalStrain_zz, "none");
        _mydata->auxUpdateDB->addValue("total_strain_xy", totalStrain_xy, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_xx", viscousStrain_xx, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_yy", viscousStrain_yy, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_zz", viscousStrain_zz, "none");
        _mydata->auxUpdateDB->addValue("viscous_strain_xy", viscousStrain_xy, "none");
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

        _mymaterial->label("Isotropic Linear Maxwell Plane Strain");
        _mymaterial->id(24);

    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain
const double pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::SMALL = 1.0e-5;

const pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::AuxConstants pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain::constants = {
    1.0e-4, // a
    2.5e-4, // b
    3.0e-4, // c
    5.0e+7, // t
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
            {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
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
            {2, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
            {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
            {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
            {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // maxwell_time
            {1, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // total_strain
            {1, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // viscous_strain
            {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // body_force
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
            {3, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
            {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
            {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
            {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // maxwell_time
            {1, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // total_strain
            {1, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // viscous_strain
            {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP3);

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
            {4, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            {0, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
            {0, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
            {0, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
            {0, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // maxwell_time
            {1, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // total_strain
            {1, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // viscous_strain
            {0, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP4
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_TriP4);

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
            {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
            {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
            {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
            {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // maxwell_time
            {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // total_strain
            {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // viscous_strain
            {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // body_force
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
            {2, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
            {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
            {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
            {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // maxwell_time
            {1, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // total_strain
            {1, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // viscous_strain
            {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // body_force
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
            {3, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
            {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
            {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
            {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // maxwell_time
            {1, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // total_strain
            {1, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // viscous_strain
            {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ3);


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
            {4, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
            {0, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
            {0, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
            {0, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
            {0, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // maxwell_time
            {1, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // total_strain
            {1, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // viscous_strain
            {0, 4, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // body_force
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ4
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_LinearStrain_QuadQ4);


// End of file
