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

#include "RigidBodyAcc2D.hh" // Implementation of test data

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    class _RigidBodyAcc2D;
} // pylith

class pylith::_RigidBodyAcc2D {
    static const double LENGTH_SCALE;
    static const double TIME_SCALE;
    static const double PRESSURE_SCALE;
    static const double TIME_SNAPSHOT; // nondimensional

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
        return 1.0;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double t) {
        return 1.5;
    } // disp_y

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    // Velocity
    static double vel_x(const double x,
                        const double y,
                        const double t) {
        return 2.2;
    } // vel_x

    static double vel_y(const double x,
                        const double y,
                        const double t) {
        return 2.7;
    } // vel_y

    static const char* vel_units(void) {
        return "m/s";
    } // vel_units

    // Acceleration
    static double acc_x(const double x,
                        const double y,
                        const double t) {
        return -2.3;
    } // vel_x

    static double acc_y(const double x,
                        const double y,
                        const double t) {
        return -2.8;
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
        assert(2 == spaceDim);
        assert(x);
        assert(2 == numComponents);
        assert(s);

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
        assert(2 == spaceDim);
        assert(x);
        assert(2 == numComponents);
        assert(s);

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
        assert(2 == spaceDim);
        assert(x);
        assert(2 == numComponents);
        assert(s);

        s[0] = acc_x(x[0], x[1], t);
        s[1] = acc_y(x[0], x[1], t);

        return 0;
    } // solnkernel_acc

    static void mmsBodyForceKernel(const PylithInt dim,
                                   const PylithInt numS,
                                   const PylithInt numA,
                                   const PylithInt sOff[],
                                   const PylithInt sOff_x[],
                                   const PylithScalar s[],
                                   const PylithScalar s_t[],
                                   const PylithScalar s_x[],
                                   const PylithInt aOff[],
                                   const PylithInt aOff_x[],
                                   const PylithScalar a[],
                                   const PylithScalar a_t[],
                                   const PylithScalar a_x[],
                                   const PylithReal t,
                                   const PylithScalar x[],
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar f0[]) {
        assert(2 == dim);
        assert(x);
        assert(f0);
        const double velocityScale = LENGTH_SCALE / TIME_SCALE;
        const double densityScale = PRESSURE_SCALE / (velocityScale * velocityScale);
        f0[0] += density(x[0], x[1]) / densityScale * acc_x(x[0], x[1], t);
        f0[1] += density(x[0], x[1]) / densityScale * acc_y(x[0], x[1], t);
    } // mmsBodyForceKernel

public:

    static
    TestLinearElasticity_Data* createData(void) {
        TestLinearElasticity_Data* data = new TestLinearElasticity_Data();assert(data);

        data->journalName = "RigidBodyAcc2D";
        data->isJacobianLinear = true;
        data->allowZeroResidual = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        data->normalizer.setLengthScale(LENGTH_SCALE);
        data->normalizer.setTimeScale(TIME_SCALE);
        data->normalizer.setPressureScale(PRESSURE_SCALE);
        data->normalizer.computeDensityScale();
        data->formulation = pylith::problems::Physics::DYNAMIC;

        // solnDiscretizations set in derived class.

        // Material information
        data->numAuxSubfields = 3;
        static const char* _auxSubfields[3] = {"density", "shear_modulus", "bulk_modulus"};
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization const*>(_auxDiscretizations);

        data->auxDB.addValue("density", density, density_units());
        data->auxDB.addValue("vp", vp, vp_units());
        data->auxDB.addValue("vs", vs, vs_units());
        data->auxDB.setCoordSys(data->cs);

        data->material.setFormulation(data->formulation);
        data->material.useBodyForce(false);
        data->rheology.useReferenceState(false);

        data->material.setDescription("Isotropic Linear Elasticity Plane Strain");
        data->material.setLabelValue(24);

        std::vector<pylith::feassemble::IntegratorDomain::ResidualKernels> mmsKernels(1);
        mmsKernels[0] = pylith::feassemble::IntegratorDomain::ResidualKernels("velocity", pylith::feassemble::Integrator::RHS, mmsBodyForceKernel, NULL);
        data->material.setMMSBodyForceKernels(mmsKernels);

        // Boundary conditions
        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        data->bcs.resize(2);
        pylith::bc::DirichletUserFn* bc = NULL;
        bc = new pylith::bc::DirichletUserFn();assert(bc);
        bc->setSubfieldName("displacement");
        bc->setLabelName("boundary");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(constrainedDOF, numConstrained);
        bc->setUserFn(solnkernel_disp);
        bc->setUserFnDot(solnkernel_vel);
        data->bcs[0] = bc;

        bc = new pylith::bc::DirichletUserFn();assert(bc);
        bc->setSubfieldName("velocity");
        bc->setLabelName("boundary");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(constrainedDOF, numConstrained);
        bc->setUserFn(solnkernel_vel);
        bc->setUserFnDot(solnkernel_acc);
        data->bcs[1] = bc;

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[2] = {
            solnkernel_disp,
            solnkernel_vel,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        static const pylith::testing::MMSTest::solution_fn _exactSolnDotFns[2] = {
            solnkernel_vel,
            solnkernel_acc,
        };
        data->exactSolnDotFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnDotFns);

        return data;
    } // createData

}; // _RigidBodyAcc2D

const double pylith::_RigidBodyAcc2D::LENGTH_SCALE = 1.0e+3;
const double pylith::_RigidBodyAcc2D::TIME_SCALE = 15.0;
const double pylith::_RigidBodyAcc2D::PRESSURE_SCALE = 2.0e+10;

// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::RigidBodyAcc2D::TriP1(void) {
    TestLinearElasticity_Data* data = pylith::_RigidBodyAcc2D::createData();assert(data);

    data->meshFilename = "data/tri.msh";
    data->useAsciiMesh = false;

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP1


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::RigidBodyAcc2D::TriP2(void) {
    TestLinearElasticity_Data* data = pylith::_RigidBodyAcc2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP2


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::RigidBodyAcc2D::TriP3(void) {
    TestLinearElasticity_Data* data = pylith::_RigidBodyAcc2D::createData();assert(data);

    data->meshFilename = "data/tri.msh";
    data->useAsciiMesh = false;

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(3, 3), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP3


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::RigidBodyAcc2D::TriP4(void) {
    TestLinearElasticity_Data* data = pylith::_RigidBodyAcc2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 4), // density
        pylith::topology::Field::Discretization(0, 4), // shear_modulus
        pylith::topology::Field::Discretization(0, 4), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(4, 4), // disp
        pylith::topology::Field::Discretization(4, 4), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP4


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::RigidBodyAcc2D::QuadQ1(void) {
    TestLinearElasticity_Data* data = pylith::_RigidBodyAcc2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ1


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::RigidBodyAcc2D::QuadQ1Distorted(void) {
    TestLinearElasticity_Data* data = pylith::_RigidBodyAcc2D::createData();assert(data);

    data->meshFilename = "data/quad_distorted.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ1Distorted


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::RigidBodyAcc2D::QuadQ2(void) {
    TestLinearElasticity_Data* data = pylith::_RigidBodyAcc2D::createData();assert(data);

    data->meshFilename = "data/quad.msh";
    data->useAsciiMesh = false;

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ2


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::RigidBodyAcc2D::QuadQ3(void) {
    TestLinearElasticity_Data* data = pylith::_RigidBodyAcc2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(3, 3), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ3


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::RigidBodyAcc2D::QuadQ4(void) {
    TestLinearElasticity_Data* data = pylith::_RigidBodyAcc2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 4), // density
        pylith::topology::Field::Discretization(0, 4), // shear_modulus
        pylith::topology::Field::Discretization(0, 4), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(4, 4), // disp
        pylith::topology::Field::Discretization(4, 4), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ4


// End of file
