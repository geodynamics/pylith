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

#include "TestIsotropicLinearElasticity.hh" // ISA TestIsotropicLinearElasticity3D

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

namespace pylith {
    namespace mmstests {
        class TestIsotropicLinearElasticity3D_UniformStrain;

        class TestIsotropicLinearElasticity3D_UniformStrain_TetP1;
        class TestIsotropicLinearElasticity3D_UniformStrain_TetP2;
        class TestIsotropicLinearElasticity3D_UniformStrain_TetP3;
        class TestIsotropicLinearElasticity3D_UniformStrain_TetP4;

        class TestIsotropicLinearElasticity3D_UniformStrain_HexQ1;
        class TestIsotropicLinearElasticity3D_UniformStrain_HexQ2;
        class TestIsotropicLinearElasticity3D_UniformStrain_HexQ3;
        class TestIsotropicLinearElasticity3D_UniformStrain_HexQ4;
    } // tests/mmstests
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain :
    public pylith::mmstests::TestIsotropicLinearElasticity {
    // Spatial database user functions for auxiiliary subfields (includes derived fields).

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

    // Solution subfields.

    static double strain_xx(void) {
        return 0.1;
    } // strain_xx

    static double strain_yy(void) {
        return 0.25;
    } // strain_yy

    static double strain_zz(void) {
        return -0.05;
    } // strain_zz

    static double strain_xy(void) {
        return 0.3;
    } // strain_xy

    static double strain_yz(void) {
        return 0;
    } // strain_yz

    static double strain_xz(void) {
        return -0.08;
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
        return strain_xy()*x + strain_yy()*y * strain_yz()*z;
    } // disp_y

    static double disp_z(const double x,
                         const double y,
                         const double z) {
        return strain_xz()*x + strain_yz()*y * strain_zz()*z;
    } // disp_z

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        CPPUNIT_ASSERT(3 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(3 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = disp_x(x[0], x[1], x[2]);
        s[1] = disp_y(x[0], x[1], x[2]);
        s[2] = disp_z(x[0], x[1], x[2]);

        return 0;
    } // solnkernel_disp

protected:

    void setUp(void) {
        TestIsotropicLinearElasticity::setUp();

        // Overwrite component names for control of journals at test level.
        GenericComponent::setName("TestIsotropicLinearElasticity3D_UniformStrain");

        CPPUNIT_ASSERT(!_data);
        _data = new TestElasticity_Data();CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->spaceDim = 3;
        _data->meshFilename = ":UNKNOWN:"; // Set in child class.
        _data->boundaryLabel = "boundary";

        CPPUNIT_ASSERT(!_data->cs);
        _data->cs = new spatialdata::geocoords::CSCart;CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->spaceDim);

        CPPUNIT_ASSERT(_data->normalizer);
        _data->normalizer->setLengthScale(1.0e+03);
        _data->normalizer->setTimeScale(2.0);
        _data->normalizer->setPressureScale(2.25e+10);
        _data->normalizer->computeDensityScale();

        _data->startTime = 0.0;
        _data->endTime = 0.1;
        _data->timeStep = 0.05;

        // solnDiscretizations set in derived class.

        _data->numAuxSubfields = 3;
        static const char* _auxSubfields[3] = {"density", "shear_modulus", "bulk_modulus"};
        _data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        CPPUNIT_ASSERT(_data->auxDB);
        _data->auxDB->addValue("density", density, density_units());
        _data->auxDB->addValue("vp", vp, vp_units());
        _data->auxDB->addValue("vs", vs, vs_units());
        _data->auxDB->setCoordSys(*_data->cs);

        CPPUNIT_ASSERT(_material);
        _material->setFormulation(pylith::problems::Physics::QUASISTATIC);
        _material->useBodyForce(false);
        _rheology->useReferenceState(false);

        _material->setDescriptiveLabel("Isotropic Linear Elascitity");
        _material->setMaterialId(24);

        static const PylithInt constrainedDOF[3] = { 0, 1, 2 };
        static const PylithInt numConstrained = 3;
        _bc->setConstrainedDOF(constrainedDOF, numConstrained);
        _bc->setMarkerLabel("boundary");
        _bc->setSubfieldName("displacement");
        _bc->setUserFn(solnkernel_disp);

    } // setUp

    // Set exact solution in domain.
    void _setExactSolution(void) {
        CPPUNIT_ASSERT(_solution);

        PetscErrorCode err = 0;
        PetscDS prob = NULL;
        err = DMGetDS(_solution->getDM(), &prob);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_disp, NULL);CPPUNIT_ASSERT(!err);
    } // _setExactSolution

}; // TestIsotropicLinearElasticity3D_UniformStrain

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_TetP1 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_TetP1,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tet.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_UniformStrain_TetP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_TetP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_TetP2 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_TetP2,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tet.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_UniformStrain_TetP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_TetP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_TetP3 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_TetP3,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tet.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_UniformStrain_TetP3
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_TetP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_TetP4 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_TetP4,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tet.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_UniformStrain_TetP4
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_TetP4);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_HexQ1 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_HexQ1,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/hex.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_UniformStrain_HexQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_HexQ1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_HexQ2 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_HexQ2,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/hex.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_UniformStrain_HexQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_HexQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_HexQ3 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_HexQ3,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/hex.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_UniformStrain_HexQ3
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_HexQ3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_HexQ4 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_UniformStrain_HexQ4,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_UniformStrain::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/hex.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_UniformStrain_HexQ4
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_UniformStrain_HexQ4);

// End of file
