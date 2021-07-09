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
        class TestIsotropicLinearElasticity2D_GravityRefState;

        class TestIsotropicLinearElasticity2D_GravityRefState_TriP1;
        class TestIsotropicLinearElasticity2D_GravityRefState_TriP2;
        class TestIsotropicLinearElasticity2D_GravityRefState_TriP3;
        class TestIsotropicLinearElasticity2D_GravityRefState_TriP4;

        class TestIsotropicLinearElasticity2D_GravityRefState_QuadQ1;
        class TestIsotropicLinearElasticity2D_GravityRefState_QuadQ2;
        class TestIsotropicLinearElasticity2D_GravityRefState_QuadQ3;
        class TestIsotropicLinearElasticity2D_GravityRefState_QuadQ4;

    } // materials
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState :
    public pylith::mmstests::TestIsotropicLinearElasticity {
    /// Spatial database user functions for auxiiliary subfields (includes derived fields).
    static const double GACC;
    static const double YMAX;

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

    static double referenceMeanStress(const double x,
                                      const double y) {
        return density(x,y) * GACC * (y-YMAX);
    } // referenceMeanStress

    static double referenceShearStress(const double x,
                                       const double y) {
        return 0.0;
    } // referenceShearStress

    static const char* stress_units(void) {
        return "Pa";
    } // stress_units

    static double referenceStrain(const double x,
                                  const double y) {
        return 0.0;
    } // referencStrain

    static const char* strain_units(void) {
        return "none";
    } // strain_units

    static double setGravityAcc_x(const double x,
                                  const double y) {
        return 0.0;
    } // setGravityAcc_x

    static double setGravityAcc_y(const double x,
                                  const double y) {
        return -GACC;
    } // setGravityAcc_y

    static const char* acc_units(void) {
        return "m/s**2";
    } // acc_units

    // Solution subfields (nondimensional).

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return 0.0;
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return 0.0;
    } // disp_y

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);
        CPPUNIT_ASSERT(x);

        s[0] = disp_x(x[0], x[1]);
        s[1] = disp_y(x[0], x[1]);

        return 0;
    } // solnkernel_disp

protected:

    void setUp(void) {
        TestIsotropicLinearElasticity::setUp();

        // Overwrite component names for control of journals at test level.
        GenericComponent::setName("TestIsotropicLinearElasticity2D_GravityRefState");
        _disableFiniteDifferenceCheck = true;

        CPPUNIT_ASSERT(!_data);
        _data = new TestElasticity_Data();CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->spaceDim = 2;
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
        _data->gravityField->setGravityDir(0.0, -1.0, 0.0);
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
        _data->auxDB->addValue("reference_strain_xx", referenceStrain, strain_units());
        _data->auxDB->addValue("reference_strain_yy", referenceStrain, strain_units());
        _data->auxDB->addValue("reference_strain_zz", referenceStrain, strain_units());
        _data->auxDB->addValue("reference_strain_xy", referenceStrain, strain_units());
        _data->auxDB->setCoordSys(*_data->cs);

        CPPUNIT_ASSERT(_material);
        _material->setFormulation(pylith::problems::Physics::QUASISTATIC);
        _material->useBodyForce(false);
        _rheology->useReferenceState(true);

        _material->setDescriptiveLabel("Isotropic Linear Elascitity Plane Strain");
        _material->setMaterialId(24);

        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
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

}; // TestIsotropicLinearElasticity2D_GravityRefState
const double pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState::GACC = 9.80665;
const double pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState::YMAX = +4.0e+3;

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState_TriP1 :
    public pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity2D_GravityRefState_TriP1,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity2D_GravityRefState::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

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

}; // TestIsotropicLinearElasticity2D_GravityRefState_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState_TriP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState_TriP2 :
    public pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity2D_GravityRefState_TriP2,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity2D_GravityRefState::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

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

}; // TestIsotropicLinearElasticity2D_GravityRefState_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState_TriP3 :
    public pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity2D_GravityRefState_TriP3,
                           TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity2D_GravityRefState::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

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

}; // TestIsotropicLinearElasticity2D_GravityRefState_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState_QuadQ1 :
    public pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity2D_GravityRefState_QuadQ1,  TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity2D_GravityRefState::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

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

}; // TestIsotropicLinearElasticity2D_GravityRefState_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState_QuadQ1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState_QuadQ2 :
    public pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity2D_GravityRefState_QuadQ2,  TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity2D_GravityRefState::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

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

}; // TestIsotropicLinearElasticity2D_GravityRefState_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState_QuadQ3 :
    public pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity2D_GravityRefState_QuadQ3,  TestIsotropicLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticity2D_GravityRefState::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

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

}; // TestIsotropicLinearElasticity2D_GravityRefState_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearElasticity2D_GravityRefState_QuadQ3);

// End of file
