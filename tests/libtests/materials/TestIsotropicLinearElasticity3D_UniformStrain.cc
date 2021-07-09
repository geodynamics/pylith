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

#include "TestIsotropicLinearElasticity3D.hh" // Implementation of cases

#include "pylith/materials/IsotropicLinearElasticity3D.hh" // USES IsotropicLinearElasticity3D
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// forward declarations
namespace pylith {
    namespace materials {
        class TestIsotropicLinearElasticity3D_UniformStrain;

        class TestIsotropicLinearElasticity3D_UniformStrain_TetP1;
        class TestIsotropicLinearElasticity3D_UniformStrain_TetP2;
        class TestIsotropicLinearElasticity3D_UniformStrain_TetP3;
        class TestIsotropicLinearElasticity3D_UniformStrain_TetP4;

        class TestIsotropicLinearElasticity3D_UniformStrain_HexQ1;
        class TestIsotropicLinearElasticity3D_UniformStrain_HexQ2;
        class TestIsotropicLinearElasticity3D_UniformStrain_HexQ3;
        class TestIsotropicLinearElasticity3D_UniformStrain_HexQ4;

    } // materials
} // pylith

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain :
    public pylith::materials::TestIsotropicLinearElasticity3D {
    /// Spatial database user functions for auxiiliary subfields (includes derived fields).
    static const double SMALL;

    // Density
    static double density(const double x,
                          const double y,
                          const double z) {
        return 2500.0;
    } // density

    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    // Vs
    static double vs(const double x,
                     const double y,
                     const double z) {
        return 3000.0;
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

    // Spatial database user functions for solution subfields.

    static double strain_xx(void) {
        return 0.1;
    } // strain_xx

    static double strain_yy(void) {
        return 0.25;
    } // strain_yy

    static double strain_zz(void) {
        return 0.15;
    } // strain_zz

    static double strain_xy(void) {
        return 0.3;
    } // strain_xy

    static double strain_yz(void) {
        return 0.35;
    } // strain_yz

    static double strain_xz(void) {
        return 0.4;
    } // strain_xz

    // Displacement
    static double disp_x(const double x,
                         const double y,
                         const double z) {
        return strain_xx()*x + strain_xy()*y + strain_xz()*z;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double z) {
        return strain_xy()*x + strain_yy()*y + strain_yz()*z;
    } // disp_y

    static double disp_z(const double x,
                         const double y,
                         const double z) {
        return strain_xz()*x + strain_yz()*y + strain_zz()*z;
    } // disp_z

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    static double disp_dot_x(const double x,
                             const double y,
                             const double z) {
        return 0.0;
    } // disp_dot_x

    static double disp_dot_y(const double x,
                             const double y,
                             const double z) {
        return 0.0;
    } // disp_dot_y

    static double disp_dot_z(const double x,
                             const double y,
                             const double z) {
        return 0.0;
    } // disp_dot_z

    static const char* disp_dot_units(void) {
        return "m/s";
    } // disp_dot_units

    // Displacement + perturbation
    static double disp_perturb_x(const double x,
                                 const double y,
                                 const double z) {
        const double perturbation = SMALL * disp_x(z, y, x);
        return disp_x(x, y, z) + perturbation;
    } // disp_perturb_x

    static double disp_perturb_y(const double x,
                                 const double y,
                                 const double z) {
        const double perturbation = SMALL * disp_y(z, y, x);
        return disp_y(x, y, z) + perturbation;
    } // disp_perturb_y

    static double disp_perturb_z(const double x,
                                 const double y,
                                 const double z) {
        const double perturbation = SMALL * disp_z(z, y, x);
        return disp_z(x, y, z) + perturbation;
    } // disp_perturb_z

protected:

    void setUp(void) {
        TestIsotropicLinearElasticity3D::setUp();
        _mydata = new TestIsotropicLinearElasticity3D_Data();CPPUNIT_ASSERT(_mydata);

        // dimension set in base class.
        // meshFilename set in derived class.
        _mydata->boundaryLabel = "boundary";

        CPPUNIT_ASSERT(_mydata->normalizer);
        _mydata->normalizer->setLengthScale(1.0e+03);
        _mydata->normalizer->setTimeScale(2.0);
        _mydata->normalizer->setDensityScale(3.0e+3);
        _mydata->normalizer->setPressureScale(2.25e+10);

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
        _material->setFormulation(pylith::problems::Physics::QUASISTATIC);
        _mymaterial->useBodyForce(false);
        _mymaterial->useReferenceState(false);

        _mymaterial->setLabel("Isotropic Linear Elascitity 3D");
        _mymaterial->id(24);
    } // setUp

}; // TestIsotropicLinearElasticity3D_UniformStrain
const double pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain::SMALL = 0.1;

// ----------------------------------------------------------------------

class pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_TetP1 :
    public pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_TetP1,
                           TestIsotropicLinearElasticity3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tet_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticity3D_UniformStrain_TetP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_TetP1);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_TetP2 :
    public pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_TetP2,
                           TestIsotropicLinearElasticity3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tet_small.mesh";

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

}; // TestIsotropicLinearElasticity3D_UniformStrain_TetP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_TetP2);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_TetP3 :
    public pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_TetP3,
                           TestIsotropicLinearElasticity3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tet_small.mesh";

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

}; // TestIsotropicLinearElasticity3D_UniformStrain_TetP3
// Leave this one out for now since it takes too long.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_TetP3);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_TetP4 :
    public pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_TetP4,
                           TestIsotropicLinearElasticity3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tet_small.mesh";

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

}; // TestIsotropicLinearElasticity3D_UniformStrain_TetP4
// Leave this one out for now since it takes too long.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_TetP4);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_HexQ1 :
    public pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_HexQ1,
                           TestIsotropicLinearElasticity3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/hex_small.mesh";

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

}; // TestIsotropicLinearElasticity3D_UniformStrain_HexQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_HexQ1);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_HexQ2 :
    public pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_HexQ2,
                           TestIsotropicLinearElasticity3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/hex_small.mesh";

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

}; // TestIsotropicLinearElasticity3D_UniformStrain_HexQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_HexQ2);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_HexQ3 :
    public pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_HexQ3,
                           TestIsotropicLinearElasticity3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/hex_small.mesh";

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

}; // TestIsotropicLinearElasticity3D_UniformStrain_HexQ3
// Leave this one out for now since it takes too long.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_HexQ3);

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_HexQ4 :
    public pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_HexQ4,
                           TestIsotropicLinearElasticity3D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/hex_small.mesh";

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

}; // TestIsotropicLinearElasticity3D_UniformStrain_HexQ4
// Leave this one out for now since it takes too long.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticity3D_UniformStrain_HexQ4);

// End of file
