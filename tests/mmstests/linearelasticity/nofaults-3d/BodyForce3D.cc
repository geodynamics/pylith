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

namespace pylith {
    namespace mmstests {
        class BodyForce3D;

        class BodyForce3D_TetP2;
        class BodyForce3D_TetP3;
        class BodyForce3D_TetP4;

        class BodyForce3D_HexQ2;
        class BodyForce3D_HexQ3;
        class BodyForce3D_HexQ4;

    } // materials
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::BodyForce3D :
    public pylith::mmstests::TestLinearElasticity {
    static const double LENGTHSCALE;
    static const double TIMESCALE;
    static const double PRESSURESCALE;
    static const double BODYFORCE;
    static const double XMAX;

    /// Spatial database user functions for auxiiliary subfields (includes derived fields).

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

    static double bodyforce_x(const double x,
                              const double y,
                              const double z) {
        return BODYFORCE;
    } // bodyforce_x

    static double bodyforce_y(const double x,
                              const double y,
                              const double z) {
        return 0.0;
    } // bodyforce_y

    static double bodyforce_z(const double x,
                              const double y,
                              const double z) {
        return 0.0;
    } // bodyforce_z

    static const char* bodyforce_units(void) {
        return "kg/(m**2*s**2)";
    } // bodyforce_units

    // Solution subfields (nondimensional)

    // Displacement
    static double disp_x(const double x,
                         const double y,
                         const double z) {
        const double velocityScale = LENGTHSCALE / TIMESCALE;
        const double densityScale = PRESSURESCALE / (velocityScale * velocityScale);
        const double accelerationScale = LENGTHSCALE / (TIMESCALE * TIMESCALE);
        const double forceScale = densityScale * accelerationScale;
        const double bodyforceN = BODYFORCE / forceScale;
        const double muN = density(x,y,z) * vs(x,y,z) * vs(x,y,z) / PRESSURESCALE;
        const double lambdaN = density(x,y,z) * vp(x,y,z) * vp(x,y,z) / PRESSURESCALE - 2.0*muN;
        const double xp = x - XMAX / LENGTHSCALE;
        return -0.5 * bodyforceN / (lambdaN + 2.0*muN) * (xp*xp);
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

        s[0] = disp_x(x[0], x[1], x[2]);
        s[1] = disp_y(x[0], x[1], x[2]);
        s[2] = disp_z(x[0], x[1], x[2]);

        return 0;
    } // solnkernel_disp

protected:

    void setUp(void) {
        TestLinearElasticity::setUp();

        // Overwrite component names for control of journals at test level.
        GenericComponent::setName("BodyForce3D");

        CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->meshFilename = ":UNKNOWN:"; // Set in child class.
        _data->boundaryLabel = "boundary";

        _data->normalizer.setLengthScale(LENGTHSCALE);
        _data->normalizer.setTimeScale(TIMESCALE);
        _data->normalizer.setPressureScale(PRESSURESCALE);
        _data->normalizer.computeDensityScale();

        // solnDiscretizations set in derived class.

        // Material information
        _data->numAuxSubfields = 4;
        static const char* _auxSubfields[4] = { // order must match order of subfields in auxiliary field
            "density",
            "body_force",
            "shear_modulus",
            "bulk_modulus",
        };
        _data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // body_force
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->auxDB.addValue("density", density, density_units());
        _data->auxDB.addValue("vp", vp, vp_units());
        _data->auxDB.addValue("vs", vs, vs_units());
        _data->auxDB.addValue("body_force_x", bodyforce_x, bodyforce_units());
        _data->auxDB.addValue("body_force_y", bodyforce_y, bodyforce_units());
        _data->auxDB.addValue("body_force_z", bodyforce_z, bodyforce_units());
        _data->auxDB.setCoordSys(_data->cs);

        _data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
        _data->material.useBodyForce(true);
        _data->rheology.useReferenceState(false);

        _data->material.setDescription("Isotropic Linear Elascitity Plane Strain");
        _data->material.setLabelValue(24);

        static const PylithInt constrainedDOF[3] = { 0, 1, 2 };
        static const PylithInt numConstrained = 3;
        _data->bc.setSubfieldName("displacement");
        _data->bc.setLabelName("boundary");
        _data->bc.setLabelValue(1);
        _data->bc.setConstrainedDOF(constrainedDOF, numConstrained);
        _data->bc.setUserFn(solnkernel_disp);

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

}; // BodyForce3D
const double pylith::mmstests::BodyForce3D::LENGTHSCALE = 1.0e+3;
const double pylith::mmstests::BodyForce3D::TIMESCALE = 2.0;
const double pylith::mmstests::BodyForce3D::PRESSURESCALE = 2.25e+10;
const double pylith::mmstests::BodyForce3D::BODYFORCE = 5.0e+3;
const double pylith::mmstests::BodyForce3D::XMAX = 4.0e+3;

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::BodyForce3D_TetP2 :
    public pylith::mmstests::BodyForce3D {
    CPPUNIT_TEST_SUB_SUITE(BodyForce3D_TetP2,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        BodyForce3D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tet.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // body_force
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // BodyForce3D_TetP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::BodyForce3D_TetP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::BodyForce3D_TetP3 :
    public pylith::mmstests::BodyForce3D {
    CPPUNIT_TEST_SUB_SUITE(BodyForce3D_TetP3,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        BodyForce3D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tet.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // body_force
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // BodyForce3D_TetP3
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::BodyForce3D_TetP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::BodyForce3D_HexQ2 :
    public pylith::mmstests::BodyForce3D {
    CPPUNIT_TEST_SUB_SUITE(BodyForce3D_HexQ2,  TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        BodyForce3D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/hex.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // body_force
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // BodyForce3D_HexQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::BodyForce3D_HexQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::BodyForce3D_HexQ3 :
    public pylith::mmstests::BodyForce3D {
    CPPUNIT_TEST_SUB_SUITE(BodyForce3D_HexQ3,  TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        BodyForce3D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/hex.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // body_force
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // BodyForce3D_HexQ3
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::BodyForce3D_HexQ3);

// End of file
