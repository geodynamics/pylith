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

#include "TestIsotropicLinearElasticity.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

namespace pylith {
    namespace mmstests {
        class TestIsotropicLinearElasticity3D_GravityRefState;

        class TestIsotropicLinearElasticity3D_GravityRefState_TetP1;
        class TestIsotropicLinearElasticity3D_GravityRefState_TetP2;
        class TestIsotropicLinearElasticity3D_GravityRefState_TetP3;
        class TestIsotropicLinearElasticity3D_GravityRefState_TetP4;

        class TestIsotropicLinearElasticity3D_GravityRefState_HexQ1;
        class TestIsotropicLinearElasticity3D_GravityRefState_HexQ2;
        class TestIsotropicLinearElasticity3D_GravityRefState_HexQ3;
        class TestIsotropicLinearElasticity3D_GravityRefState_HexQ4;

    } // materials
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState :
    public pylith::mmstests::TestIsotropicLinearElasticity {
    /// Spatial database user functions for auxiiliary subfields (includes derived fields).
    static const double GACC;
    static const double ZMAX;

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

    static double referenceMeanStress(const double x,
                                      const double y,
                                      const double z) {
        return density(x,y,z) * GACC * (z-ZMAX);
    } // referenceMeanStress

    static double referenceShearStress(const double x,
                                       const double y,
                                       const double z) {
        return 0.0;
    } // referenceShearStress

    static const char* stress_units(void) {
        return "Pa";
    } // stress_units

    static double referenceStrain(const double x,
                                  const double y,
                                  const double z) {
        return 0.0;
    } // referencStrain

    static const char* strain_units(void) {
        return "none";
    } // strain_units

    // Solution subfields (nondimensional).

    // Displacement
    static double disp_x(const double x,
                         const double y,
                         const double z) {
        return 0.0;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double z) {
        return 0.0;
    } // disp_y

    static double disp_z(const double x,
                         const double y,
                         const double z) {
        return 0.0;
    } // disp_z

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        CPPUNIT_ASSERT(3 == spaceDim);
        CPPUNIT_ASSERT(3 == numComponents);
        CPPUNIT_ASSERT(s);
        CPPUNIT_ASSERT(x);

        s[0] = disp_x(x[0], x[1], x[2]);
        s[1] = disp_y(x[0], x[1], x[2]);
        s[2] = disp_y(x[0], x[1], x[2]);

        return 0;
    } // solnkernel_disp

protected:

    void setUp(void) {
        TestIsotropicLinearElasticity::setUp();

        // Overwrite component names for control of journals at test level.
        GenericComponent::setName("TestIsotropicLinearElasticity3D_GravityRefState");
        _disableFiniteDifferenceCheck = true;

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

        delete _data->gravityField;_data->gravityField = new spatialdata::spatialdb::GravityField();
        _data->gravityField->setGravityDir(0.0, 0.0, -1.0);
        _data->gravityField->setGravityAcc(GACC);

        _data->startTime = 0.0;
        _data->endTime = 0.1;
        _data->timeStep = 0.05;

        // solnDiscretizations set in derived class.

        _data->numAuxSubfields = 6;
        static const char* _auxSubfields[6] = { // order must match order of subfields in auxiliary field
            "density",
            "gravitational_acceleration",
            "reference_stress",
            "reference_strain",
            "shear_modulus",
            "bulk_modulus",
        };
        _data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 1), // reference_stress
            pylith::topology::Field::Discretization(0, 1), // reference_strain
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        CPPUNIT_ASSERT(_data->auxDB);
        _data->auxDB->addValue("density", density, density_units());
        _data->auxDB->addValue("vp", vp, vp_units());
        _data->auxDB->addValue("vs", vs, vs_units());
        _data->auxDB->addValue("reference_stress_xx", referenceMeanStress, stress_units());
        _data->auxDB->addValue("reference_stress_yy", referenceMeanStress, stress_units());
        _data->auxDB->addValue("reference_stress_zz", referenceMeanStress, stress_units());
        _data->auxDB->addValue("reference_stress_xy", referenceShearStress, stress_units());
        _data->auxDB->addValue("reference_stress_yz", referenceShearStress, stress_units());
        _data->auxDB->addValue("reference_stress_xz", referenceShearStress, stress_units());
        _data->auxDB->addValue("reference_strain_xx", referenceStrain, strain_units());
        _data->auxDB->addValue("reference_strain_yy", referenceStrain, strain_units());
        _data->auxDB->addValue("reference_strain_zz", referenceStrain, strain_units());
        _data->auxDB->addValue("reference_strain_xy", referenceStrain, strain_units());
        _data->auxDB->addValue("reference_strain_yz", referenceStrain, strain_units());
        _data->auxDB->addValue("reference_strain_xz", referenceStrain, strain_units());
        _data->auxDB->setCoordSys(*_data->cs);

        CPPUNIT_ASSERT(_material);
        _material->setFormulation(pylith::problems::Physics::QUASISTATIC);
        _material->useBodyForce(false);
        _rheology->useReferenceState(true);

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

}; // TestIsotropicLinearElasticity3D_GravityRefState
const double pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState::GACC = 9.80665;
const double pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState::ZMAX = +4.0e+3;

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState_TetP1 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_GravityRefState_TetP1,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_GravityRefState::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tet.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 1), // reference_stress
            pylith::topology::Field::Discretization(0, 1), // reference_strain
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_GravityRefState_TetP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState_TetP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState_TetP2 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_GravityRefState_TetP2,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_GravityRefState::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tet.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 2), // reference_stress
            pylith::topology::Field::Discretization(0, 2), // reference_strain
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_GravityRefState_TetP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState_TetP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState_TetP3 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_GravityRefState_TetP3,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_GravityRefState::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tet.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 3), // reference_stress
            pylith::topology::Field::Discretization(0, 3), // reference_strain
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_GravityRefState_TetP3
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState_TetP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState_HexQ1 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_GravityRefState_HexQ1,  TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_GravityRefState::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/hex.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 1), // reference_stress
            pylith::topology::Field::Discretization(0, 1), // reference_strain
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_GravityRefState_HexQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState_HexQ1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState_HexQ2 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_GravityRefState_HexQ2,  TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_GravityRefState::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/hex.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 2), // reference_stress
            pylith::topology::Field::Discretization(0, 2), // reference_strain
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_GravityRefState_HexQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState_HexQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState_HexQ3 :
    public pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D_GravityRefState_HexQ3,  TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity3D_GravityRefState::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/hex.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 3), // reference_stress
            pylith::topology::Field::Discretization(0, 3), // reference_strain
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearElasticity3D_GravityRefState_HexQ3
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity3D_GravityRefState_HexQ3);

// End of file
