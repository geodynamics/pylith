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

#include "TestLinearElasticity.hh" // ISA TestLinearElasticity2D

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

namespace pylith {
    namespace mmstests {
        class PlanePWave2D;

        class PlanePWave2D_TriP1;
        class PlanePWave2D_TriP2;
        class PlanePWave2D_TriP3;
        class PlanePWave2D_TriP4;

        class PlanePWave2D_QuadQ1;
        class PlanePWave2D_QuadQ1Distorted;
        class PlanePWave2D_QuadQ2;
        class PlanePWave2D_QuadQ3;
        class PlanePWave2D_QuadQ4;
    } // tests/mmstests
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::PlanePWave2D : public pylith::mmstests::TestLinearElasticity {
    // Spatial database user functions for auxiiliary subfields (includes derived fields).

    static const double LENGTH_SCALE;
    static const double TIME_SCALE;
    static const double WAVELENGTH; // (in terms of LENGTH_SCALE)
    static const double PRESSURE_SCALE;
    static const double TIME_SNAPSHOT; // nondimensional
    static const double AMPLITUDE;

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

    // Solution subfields.

    // Displacement
    static double disp_x(const double x,
                         const double y,
                         const double t) {
        const double pi = M_PI;
        const double l = WAVELENGTH;
        const double velocityScale = LENGTH_SCALE / TIME_SCALE;
        const double c = vp(x,y) / velocityScale;
        return AMPLITUDE*sin(2.0*pi*(x-c*t)/l);
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double t) {
        return 0.0;
    } // disp_y

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    // Velocity
    static double vel_x(const double x,
                        const double y,
                        const double t) {
        const double pi = M_PI;
        const double l = WAVELENGTH;
        const double velocityScale = LENGTH_SCALE / TIME_SCALE;
        const double c = vp(x,y) / velocityScale;
        return -AMPLITUDE*2.0*pi*c/l * cos(2.0*pi*(x-c*t)/l);
    } // vel_x

    static double vel_y(const double x,
                        const double y,
                        const double t) {
        return 0.0;
    } // vel_y

    static const char* vel_units(void) {
        return "m/s";
    } // vel_units

    // Acceleration
    static double acc_x(const double x,
                        const double y,
                        const double t) {
        const double pi = M_PI;
        const double l = WAVELENGTH;
        const double velocityScale = LENGTH_SCALE / TIME_SCALE;
        const double c = vp(x,y) / velocityScale;
        return -AMPLITUDE*pow(2.0*pi*c/l, 2) * sin(2.0*pi*(x-c*t)/l);
    } // vel_x

    static double acc_y(const double x,
                        const double y,
                        const double t) {
        return 0.0;
    } // vel_y

    static const char* acc_units(void) {
        return "m/s*2";
    } // vel_units

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = disp_x(x[0], x[1], t);
        s[1] = disp_y(x[0], x[1], t);

        return 0;
    } // solnkernel_disp

    static PetscErrorCode solnkernel_vel(PetscInt spaceDim,
                                         PetscReal t,
                                         const PetscReal x[],
                                         PetscInt numComponents,
                                         PetscScalar* s,
                                         void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = vel_x(x[0], x[1], t);
        s[1] = vel_y(x[0], x[1], t);

        return 0;
    } // solnkernel_vel

    static PetscErrorCode solnkernel_acc(PetscInt spaceDim,
                                         PetscReal t,
                                         const PetscReal x[],
                                         PetscInt numComponents,
                                         PetscScalar* s,
                                         void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = acc_x(x[0], x[1], t);
        s[1] = acc_y(x[0], x[1], t);

        return 0;
    } // solnkernel_acc

protected:

    void setUp(void) {
        TestLinearElasticity::setUp();

        // Overwrite component names for control of journals at test level.
        GenericComponent::setName("PlanePWave2D");

        CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->meshFilename = ":UNKNOWN:"; // Set in child class.
        _data->boundaryLabel = "boundary";

        _data->normalizer.setLengthScale(LENGTH_SCALE);
        _data->normalizer.setTimeScale(TIME_SCALE);
        _data->normalizer.setPressureScale(PRESSURE_SCALE);
        _data->normalizer.computeDensityScale();
        _data->formulation = pylith::problems::Physics::DYNAMIC;

        _data->t = TIME_SNAPSHOT;
        _data->dt = 0.05;

        // solnDiscretizations set in derived class.

        // Material information
        _data->numAuxSubfields = 3;
        static const char* _auxSubfields[3] = {"density", "shear_modulus", "bulk_modulus"};
        _data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization const*>(_auxDiscretizations);

        _data->auxDB.addValue("density", density, density_units());
        _data->auxDB.addValue("vp", vp, vp_units());
        _data->auxDB.addValue("vs", vs, vs_units());
        _data->auxDB.setCoordSys(_data->cs);

        _data->material.setFormulation(_data->formulation);
        _data->material.useBodyForce(false);
        _data->rheology.useReferenceState(false);

        _data->material.setDescription("Isotropic Linear Elasticity Plane Strain");
        _data->material.setLabelValue(24);

        // Boundary conditions
        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        pylith::bc::DirichletUserFn* bc = NULL;
        _data->bcs.resize(2);
        bc = new pylith::bc::DirichletUserFn();CPPUNIT_ASSERT(bc);
        bc->setSubfieldName("displacement");
        bc->setLabelName("boundary");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(constrainedDOF, numConstrained);
        bc->setUserFn(solnkernel_disp);
        bc->setUserFnDot(solnkernel_vel);
        _data->bcs[0] = bc;

        bc = new pylith::bc::DirichletUserFn();CPPUNIT_ASSERT(bc);
        bc->setSubfieldName("velocity");
        bc->setLabelName("boundary");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(constrainedDOF, numConstrained);
        bc->setUserFn(solnkernel_vel);
        bc->setUserFnDot(solnkernel_acc);
        _data->bcs[1] = bc;
    } // setUp

    // Set exact solution in domain.
    void _setExactSolution(void) {
        const pylith::topology::Field* solution = _problem->getSolution();
        CPPUNIT_ASSERT(solution);

        PetscErrorCode err = 0;
        PetscDS ds = NULL;
        err = DMGetDS(solution->getDM(), &ds);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(ds, 0, solnkernel_disp, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(ds, 1, solnkernel_vel, NULL);CPPUNIT_ASSERT(!err);

        err = PetscDSSetExactSolutionTimeDerivative(ds, 0, solnkernel_vel, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolutionTimeDerivative(ds, 1, solnkernel_acc, NULL);CPPUNIT_ASSERT(!err);
    } // _setExactSolution

}; // PlanePWave2D

const double pylith::mmstests::PlanePWave2D::LENGTH_SCALE = 1.0e+3;
const double pylith::mmstests::PlanePWave2D::TIME_SCALE = 10.0;
const double pylith::mmstests::PlanePWave2D::PRESSURE_SCALE = 3.0e+10;
const double pylith::mmstests::PlanePWave2D::WAVELENGTH = 1.0e+4;
const double pylith::mmstests::PlanePWave2D::TIME_SNAPSHOT = 7.657345769747113;
const double pylith::mmstests::PlanePWave2D::AMPLITUDE = 1.0e+2;

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::PlanePWave2D_TriP1 :
    public pylith::mmstests::PlanePWave2D {
    CPPUNIT_TEST_SUB_SUITE(PlanePWave2D_TriP1,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        PlanePWave2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.msh";
        _data->useAsciiMesh = false;
        _tolerance = 1.0e-4;

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(1, 1), // disp
            pylith::topology::Field::Discretization(1, 1), // vel
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // PlanePWave2D_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::PlanePWave2D_TriP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::PlanePWave2D_TriP2 :
    public pylith::mmstests::PlanePWave2D {
    CPPUNIT_TEST_SUB_SUITE(PlanePWave2D_TriP2,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        PlanePWave2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";
        _tolerance = 5.0e-7;

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(2, 2), // vel
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // PlanePWave2D_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::PlanePWave2D_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::PlanePWave2D_TriP3 :
    public pylith::mmstests::PlanePWave2D {
    CPPUNIT_TEST_SUB_SUITE(PlanePWave2D_TriP3,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        PlanePWave2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.msh";
        _data->useAsciiMesh = false;
        _tolerance = 1.0e-9;

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(3, 3), // disp
            pylith::topology::Field::Discretization(3, 3), // vel
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // PlanePWave2D_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::PlanePWave2D_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::PlanePWave2D_TriP4 :
    public pylith::mmstests::PlanePWave2D {
    CPPUNIT_TEST_SUB_SUITE(PlanePWave2D_TriP4,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        PlanePWave2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";
        _tolerance = 1.0e-9;

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(4, 4), // disp
            pylith::topology::Field::Discretization(4, 4), // vel
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // PlanePWave2D_TriP4
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::PlanePWave2D_TriP4);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::PlanePWave2D_QuadQ1 :
    public pylith::mmstests::PlanePWave2D {
    CPPUNIT_TEST_SUB_SUITE(PlanePWave2D_QuadQ1,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        PlanePWave2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";
        _tolerance = 1.0e-8;

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(1, 1), // disp
            pylith::topology::Field::Discretization(1, 1), // vel
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // PlanePWave2D_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::PlanePWave2D_QuadQ1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::PlanePWave2D_QuadQ1Distorted :
    public pylith::mmstests::PlanePWave2D {
    CPPUNIT_TEST_SUB_SUITE(PlanePWave2D_QuadQ1Distorted,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        PlanePWave2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad_distorted.mesh";
        _tolerance = 1.0e-5;

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(1, 1), // disp
            pylith::topology::Field::Discretization(1, 1), // vel
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // PlanePWave2D_QuadQ1Distorted
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::PlanePWave2D_QuadQ1Distorted);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::PlanePWave2D_QuadQ2 :
    public pylith::mmstests::PlanePWave2D {
    CPPUNIT_TEST_SUB_SUITE(PlanePWave2D_QuadQ2,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        PlanePWave2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.msh";
        _data->useAsciiMesh = false;
        _tolerance = 1.0e-9;

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(2, 2), // vel
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // PlanePWave2D_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::PlanePWave2D_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::PlanePWave2D_QuadQ3 :
    public pylith::mmstests::PlanePWave2D {
    CPPUNIT_TEST_SUB_SUITE(PlanePWave2D_QuadQ3,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        PlanePWave2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";
        _tolerance = 1.0e-9;

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(3, 3), // disp
            pylith::topology::Field::Discretization(3, 3), // vel
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // PlanePWave2D_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::PlanePWave2D_QuadQ3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::PlanePWave2D_QuadQ4 :
    public pylith::mmstests::PlanePWave2D {
    CPPUNIT_TEST_SUB_SUITE(PlanePWave2D_QuadQ4,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        PlanePWave2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(4, 4), // disp
            pylith::topology::Field::Discretization(4, 4), // vel
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // PlanePWave2D_QuadQ4
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::PlanePWave2D_QuadQ4);

// End of file
