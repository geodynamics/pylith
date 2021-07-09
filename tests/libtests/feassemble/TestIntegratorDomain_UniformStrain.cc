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

#include "TestIntegratorDomain.hh" // Implementation of cases

#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// forward declarations
namespace pylith {
    namespace materials {
        class TestIntegratorDomain_UniformStrain;

        class TestIntegratorDomain_UniformStrain_TriP1;
        class TestIntegratorDomain_UniformStrain_TriP2;

        class TestIntegratorDomain_UniformStrain_QuadQ1;
        class TestIntegratorDomain_UniformStrain_QuadQ2;

    } // materials
} // pylith

// ----------------------------------------------------------------------
class pylith::feassemble::TestIntegratorDomain_UniformStrain2D :
    public pylith::feassemble::TestIntegratorDomain {
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
        return -0.25;
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
        TestIntegratorDomain::setUp();
        _data = new TestIntegratorDomain_Data();CPPUNIT_ASSERT(_data);

        _data->dimension = 2;
        // meshFilename set in derived class.
        _data->boundaryLabel = "boundary";
        _data->materialId = 24;

        _data->cs = new spatialdata::spatialdb::CSCart();CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->dimension);

        CPPUNIT_ASSERT(_data->normalizer);
        _data->normalizer->setLengthScale(1.0e+03);
        _data->normalizer->setTimeScale(2.0);
        _data->normalizer->setDensityScale(3.0e+3);
        _data->normalizer->setPressureScale(2.25e+10);

        _data->t = 1.0;
        _data->dt = 0.05;
        _data->tindex = 20;
        _data->s_tshift = 1.0 / _data->dt;

        _data->numSolutionSubfields = 1;
        // solnDiscretizations set in derived class.

        CPPUNIT_ASSERT(_data->solutionDB);
        _data->solutionDB->addValue("displacement_x", disp_x, disp_units());
        _data->solutionDB->addValue("displacement_y", disp_y, disp_units());
        _data->solutionDB->addValue("displacement_dot_x", disp_dot_x, disp_dot_units());
        _data->solutionDB->addValue("displacement_dot_y", disp_dot_y, disp_dot_units());

        CPPUNIT_ASSERT(_data->perturbationDB);
        _data->perturbationDB->addValue("displacement_x", disp_perturb_x, disp_units());
        _data->perturbationDB->addValue("displacement_y", disp_perturb_y, disp_units());
        _data->perturbationDB->addValue("displacement_dot_x", disp_dot_x, disp_dot_units());
        _data->perturbationDB->addValue("displacement_dot_y", disp_dot_y, disp_dot_units());

        _data->numAuxiliarySubfields = 3;
        static const char* _auxiliarySubfields[3] = {"density", "shear_modulus", "bulk_modulus"};
        _data->auxiliarySubfields = _auxiliarySubfields;
        static const pylith::topology::Field::Discretization _auxiliaryDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->auxiliaryDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxiliaryDiscretizations);

        CPPUNIT_ASSERT(_data->auxiliaryDB);
        _data->auxiliaryDB->addValue("density", density, density_units());
        _data->auxiliaryDB->addValue("vp", vp, vp_units());
        _data->auxiliaryDB->addValue("vs", vs, vs_units());
        _data->auxiliaryDB->addValue("shear_modulus", shearModulus, shearModulus_units());
        _data->auxiliaryDB->addValue("bulk_modulus", bulkModulus, bulkModulus_units());

        CPPUNIT_ASSERT(_integrator);
        _integrator->useInertia(false);
        _integrator->useBodyForce(false);
        _integrator->useReferenceState(false);

        _integrator->setLabel("Isotropic Linear Elascitity Plane Strain");
        _integrator->id(24);
    } // setUp

}; // TestIntegratorDomain_UniformStrain
const double pylith::feassemble::TestIntegratorDomain_UniformStrain::SMALL = 0.1;

// ----------------------------------------------------------------------

class pylith::feassemble::TestIntegratorDomain_UniformStrain_TriP1 :
    public pylith::feassemble::TestIntegratorDomain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIntegratorDomain_UniformStrain_TriP1,
                           TestIntegratorDomain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIntegratorDomain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri_small.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIntegratorDomain_UniformStrain_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::feassemble::TestIntegratorDomain_UniformStrain_TriP1);

#if 0
// ----------------------------------------------------------------------
class pylith::feassemble::TestIntegratorDomain_UniformStrain_TriP2 :
    public pylith::feassemble::TestIntegratorDomain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIntegratorDomain_UniformStrain_TriP2,
                           TestIntegratorDomain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIntegratorDomain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri_small.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxiliaryDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxiliaryDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxiliaryDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIntegratorDomain_UniformStrain_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::feassemble::TestIntegratorDomain_UniformStrain_TriP2);

// ----------------------------------------------------------------------
class pylith::feassemble::TestIntegratorDomain_UniformStrain_TriP3 :
    public pylith::feassemble::TestIntegratorDomain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIntegratorDomain_UniformStrain_TriP3,
                           TestIntegratorDomain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIntegratorDomain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri_small.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxiliaryDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxiliaryDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxiliaryDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIntegratorDomain_UniformStrain_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::feassemble::TestIntegratorDomain_UniformStrain_TriP3);

// ----------------------------------------------------------------------
class pylith::feassemble::TestIntegratorDomain_UniformStrain_TriP4 :
    public pylith::feassemble::TestIntegratorDomain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIntegratorDomain_UniformStrain_TriP4,
                           TestIntegratorDomain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIntegratorDomain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri_small.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxiliaryDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _data->auxiliaryDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxiliaryDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIntegratorDomain_UniformStrain_TriP4
// Leave this out for now to shorten runtime.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::feassemble::TestIntegratorDomain_UniformStrain_TriP4);

// ----------------------------------------------------------------------
class pylith::feassemble::TestIntegratorDomain_UniformStrain_QuadQ1 :
    public pylith::feassemble::TestIntegratorDomain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIntegratorDomain_UniformStrain_QuadQ1,
                           TestIntegratorDomain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIntegratorDomain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad_aligned.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxiliaryDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->auxiliaryDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxiliaryDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIntegratorDomain_UniformStrain_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::feassemble::TestIntegratorDomain_UniformStrain_QuadQ1);

// ----------------------------------------------------------------------
class pylith::feassemble::TestIntegratorDomain_UniformStrain_QuadQ2 :
    public pylith::feassemble::TestIntegratorDomain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIntegratorDomain_UniformStrain_QuadQ2,
                           TestIntegratorDomain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIntegratorDomain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad_aligned.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxiliaryDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxiliaryDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxiliaryDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIntegratorDomain_UniformStrain_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::feassemble::TestIntegratorDomain_UniformStrain_QuadQ2);

// ----------------------------------------------------------------------
class pylith::feassemble::TestIntegratorDomain_UniformStrain_QuadQ3 :
    public pylith::feassemble::TestIntegratorDomain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIntegratorDomain_UniformStrain_QuadQ3,
                           TestIntegratorDomain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIntegratorDomain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad_aligned.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxiliaryDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxiliaryDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxiliaryDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIntegratorDomain_UniformStrain_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::feassemble::TestIntegratorDomain_UniformStrain_QuadQ3);

// ----------------------------------------------------------------------
class pylith::feassemble::TestIntegratorDomain_UniformStrain_QuadQ4 :
    public pylith::feassemble::TestIntegratorDomain_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIntegratorDomain_UniformStrain_QuadQ4,
                           TestIntegratorDomain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIntegratorDomain_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad_aligned.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxiliaryDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _data->auxiliaryDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxiliaryDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIntegratorDomain_UniformStrain_QuadQ4
// Leave this out for now to shorten runtime.
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::feassemble::TestIntegratorDomain_UniformStrain_QuadQ4);
#endif

// End of file
