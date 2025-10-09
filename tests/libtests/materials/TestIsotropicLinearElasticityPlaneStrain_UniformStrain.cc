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

#include "TestIsotropicLinearElasticityPlaneStrain.hh" // Implementation of cases

#include "pylith/materials/IsotropicLinearElasticityPlaneStrain.hh" // USES IsotropicLinearElasticityPlaneStrain
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Scales.hh" // USES Scales

// forward declarations
namespace pylith {
    namespace materials {
        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain;

        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP1;
        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP2;
        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP3;
        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP4;

        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ1;
        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ2;
        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ3;
        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ4;

    } // materials
} // pylith

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain {
    /// Spatial database user functions for auxiiliary subfields (includes derived fields).
    static const double SMALL;

    // Density
    static double density(const double x,
                          const double y) {
        return 2500.0;
    } // density

    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    // Vs
    static double vs(const double x,
                     const double y) {
        return 3000.0;
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

    // Spatial database user functions for solution subfields.

    static double strain_xx(void) {
        return 0.1;
    } // strain_xx

    static double strain_yy(void) {
        return 0.25;
    } // strain_yy

    static double strain_xy(void) {
        return 0.3;
    } // strain_xy

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return strain_xx()*x + strain_xy()*y;
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return strain_xy()*x + strain_yy()*y;
    } // disp_y

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    static double disp_dot_x(const double x,
                             const double y) {
        return 0.0;
    } // disp_dot_x

    static double disp_dot_y(const double x,
                             const double y) {
        return 0.0;
    } // disp_dot_y

    static const char* disp_dot_units(void) {
        return "m/s";
    } // disp_dot_units

    // Displacement + perturbation
    static double disp_perturb_x(const double x,
                                 const double y) {
        const double perturbation = SMALL * disp_x(y, x);
        return disp_x(x, y) + perturbation;
    } // disp_perturb_x

    static double disp_perturb_y(const double x,
                                 const double y) {
        const double perturbation = SMALL * disp_y(y, x);
        return disp_y(x, y) + perturbation;
    } // disp_perturb_y

protected:

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain::setUp();
        _mydata = new TestIsotropicLinearElasticityPlaneStrain_Data();CPPUNIT_ASSERT(_mydata);

        // dimension set in base class.
        // meshFilename set in derived class.
        _mydata->boundaryLabel = "boundary";

        CPPUNIT_ASSERT(_mydata->scales);
        _mydata->scales->setLengthScale(1.0);
        _mydata->scales->setTimeScale(2.0);
        _mydata->scales->setRigidityScale(2.5e+6);

        _mydata->t = 1.0;
        _mydata->dt = 0.05;
        _mydata->s_tshift = 1.0 / _mydata->dt;

        // solnDiscretizations set in derived class.

        _mydata->numAuxSubfields = 3;
        static const char* _auxSubfields[3] = {"density", "shear_modulus", "bulk_modulus"};
        _mydata->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        CPPUNIT_ASSERT(_mydata->auxDB);
        _mydata->auxDB->addValue("density", density, density_units());
        _mydata->auxDB->addValue("vp", vp, vp_units());
        _mydata->auxDB->addValue("vs", vs, vs_units());
        _mydata->auxDB->addValue("shear_modulus", shearModulus, shearModulus_units());
        _mydata->auxDB->addValue("bulk_modulus", bulkModulus, bulkModulus_units());

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
        _mymaterial->useBodyForce(false);
        _mymaterial->useReferenceState(false);

        _mymaterial->setLabel("Isotropic Linear Elascitity Plane Strain");
        _mymaterial->id(24);
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain
const double pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain::SMALL = 0.1;

// ----------------------------------------------------------------------

class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP1 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP1,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tri_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP1);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP2 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP2,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tri_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP2);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP3 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP3,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tri_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP3);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP4 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP4,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tri_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP4
// Leave this out for now to shorten runtime.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP4);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ1 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ1,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/quad_aligned.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ1);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ2 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ2,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/quad_aligned.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ2);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ3 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ3,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/quad_aligned.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ3);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ4 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ4,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/quad_aligned.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ4
// Leave this out for now to shorten runtime.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadQ4);

// End of file
