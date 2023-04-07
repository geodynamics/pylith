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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestLinearElasticity.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

namespace pylith {
    namespace mmstests {
        class GravityRefState2D;

        class GravityRefState2D_TriP1;
        class GravityRefState2D_TriP2;
        class GravityRefState2D_TriP3;
        class GravityRefState2D_TriP4;

        class GravityRefState2D_QuadQ1;
        class GravityRefState2D_QuadQ2;
        class GravityRefState2D_QuadQ3;
        class GravityRefState2D_QuadQ4;

    } // materials
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::GravityRefState2D :
    public pylith::mmstests::TestLinearElasticity {
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
        TestLinearElasticity::setUp();

        // Overwrite component names for control of journals at test level.
        GenericComponent::setName("GravityRefState2D");
        _disableFiniteDifferenceCheck = true; // Reference state matches, so no Jacobian is formed.

        CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->meshFilename = ":UNKNOWN:"; // Set in child class.
        _data->boundaryLabel = "boundary";

        _data->normalizer.setLengthScale(1.0e+03);
        _data->normalizer.setTimeScale(2.0);
        _data->normalizer.setPressureScale(2.25e+10);
        _data->normalizer.computeDensityScale();

        delete _data->gravityField;_data->gravityField = new spatialdata::spatialdb::GravityField();
        _data->gravityField->setGravityDir(0.0, -1.0, 0.0);
        _data->gravityField->setGravityAcc(GACC);

        // solnDiscretizations set in derived class.

        // Material information
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

        _data->auxDB.addValue("density", density, density_units());
        _data->auxDB.addValue("vp", vp, vp_units());
        _data->auxDB.addValue("vs", vs, vs_units());
        _data->auxDB.addValue("reference_stress_xx", referenceMeanStress, stress_units());
        _data->auxDB.addValue("reference_stress_yy", referenceMeanStress, stress_units());
        _data->auxDB.addValue("reference_stress_zz", referenceMeanStress, stress_units());
        _data->auxDB.addValue("reference_stress_xy", referenceShearStress, stress_units());
        _data->auxDB.addValue("reference_strain_xx", referenceStrain, strain_units());
        _data->auxDB.addValue("reference_strain_yy", referenceStrain, strain_units());
        _data->auxDB.addValue("reference_strain_zz", referenceStrain, strain_units());
        _data->auxDB.addValue("reference_strain_xy", referenceStrain, strain_units());
        _data->auxDB.setCoordSys(_data->cs);

        _data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
        _data->material.useBodyForce(false);
        _data->rheology.useReferenceState(true);

        _data->material.setDescription("Isotropic Linear Elascitity Plane Strain");
        _data->material.setLabelValue(24);

        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        _data->bcs.resize(1);
        pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();CPPUNIT_ASSERT(bc);
        bc->setSubfieldName("displacement");
        bc->setLabelName("boundary");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(constrainedDOF, numConstrained);
        bc->setUserFn(solnkernel_disp);
        _data->bcs[0] = bc;

    } // setUp

    // Set exact solution in domain.
    void _setExactSolution(void) {
        const pylith::topology::Field* solution = _problem->getSolution();
        CPPUNIT_ASSERT(solution);

        PetscErrorCode err = 0;
        PetscDS prob = NULL;
        err = DMGetDS(solution->getDM(), &prob);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_disp, NULL);CPPUNIT_ASSERT(!err);
    } // _setExactSolution

}; // GravityRefState2D
const double pylith::mmstests::GravityRefState2D::GACC = 9.80665;
const double pylith::mmstests::GravityRefState2D::YMAX = +4.0e+3;

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::GravityRefState2D_TriP1 :
    public pylith::mmstests::GravityRefState2D {
    CPPUNIT_TEST_SUB_SUITE(GravityRefState2D_TriP1,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        GravityRefState2D::setUp();
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

}; // GravityRefState2D_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::GravityRefState2D_TriP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::GravityRefState2D_TriP2 :
    public pylith::mmstests::GravityRefState2D {
    CPPUNIT_TEST_SUB_SUITE(GravityRefState2D_TriP2,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        GravityRefState2D::setUp();
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

}; // GravityRefState2D_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::GravityRefState2D_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::GravityRefState2D_TriP3 :
    public pylith::mmstests::GravityRefState2D {
    CPPUNIT_TEST_SUB_SUITE(GravityRefState2D_TriP3,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        GravityRefState2D::setUp();
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

}; // GravityRefState2D_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::GravityRefState2D_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::GravityRefState2D_QuadQ1 :
    public pylith::mmstests::GravityRefState2D {
    CPPUNIT_TEST_SUB_SUITE(GravityRefState2D_QuadQ1,  TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        GravityRefState2D::setUp();
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

}; // GravityRefState2D_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::GravityRefState2D_QuadQ1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::GravityRefState2D_QuadQ2 :
    public pylith::mmstests::GravityRefState2D {
    CPPUNIT_TEST_SUB_SUITE(GravityRefState2D_QuadQ2,  TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        GravityRefState2D::setUp();
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

}; // GravityRefState2D_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::GravityRefState2D_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::GravityRefState2D_QuadQ3 :
    public pylith::mmstests::GravityRefState2D {
    CPPUNIT_TEST_SUB_SUITE(GravityRefState2D_QuadQ3,  TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        GravityRefState2D::setUp();
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

}; // GravityRefState2D_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::GravityRefState2D_QuadQ3);

// End of file
