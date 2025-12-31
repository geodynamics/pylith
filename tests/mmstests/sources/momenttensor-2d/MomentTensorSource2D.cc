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

#include "MomentTensorSource2D.hh" // Implementation of test data

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "pylith/scales/ElasticityScales.hh" // USES ElasticityScales

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    class _MomentTensorSource2D;
} // pylith

class pylith::_MomentTensorSource2D {
    // Physical constants for the problem
    static constexpr double LENGTH_SCALE = 1.0e+4;     // 10 km
    static constexpr double VELOCITY_SCALE = 3000.0;   // 3 km/s (shear wave speed)
    static constexpr double TIME_SCALE = LENGTH_SCALE / VELOCITY_SCALE;
    static constexpr double DENSITY_SCALE = 2500.0;    // kg/m^3
    static constexpr double PRESSURE_SCALE = DENSITY_SCALE * VELOCITY_SCALE * VELOCITY_SCALE;

    static constexpr double TIME_SNAPSHOT = 0.5;       // Time for test evaluation
    static constexpr double CENTER_FREQUENCY = 1.0;    // Hz
    static constexpr double TIME_DELAY = 0.5;          // s

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

    // Moment tensor components (in 2D: Mxx, Myy, Mxy, Myx)
    // Using an isotropic source: Mxx = Myy = 1, Mxy = Myx = 0
    static double moment_xx(const double x, const double y) {
        return 1.0e+15; // Newton-meters (Pa*m^3)
    }
    static double moment_yy(const double x, const double y) {
        return 1.0e+15;
    }
    static double moment_xy(const double x, const double y) {
        return 0.0;
    }
    static double moment_yx(const double x, const double y) {
        return 0.0;
    }
    static const char* moment_units(void) {
        return "Pa*m**3";
    }

    // Time delay
    static double time_delay(const double x, const double y) {
        return TIME_DELAY;
    }
    static const char* time_delay_units(void) {
        return "s";
    }

    // Center frequency
    static double center_frequency(const double x, const double y) {
        return CENTER_FREQUENCY;
    }
    static const char* center_frequency_units(void) {
        return "1/s";
    }

    // Solution subfields.
    // For this MMS test, we use zero displacement/velocity to isolate testing of the
    // moment tensor source term contribution. This verifies that:
    // 1. The auxiliary fields are correctly set up
    // 2. The kernel correctly computes the source contribution
    // 3. The discretization converges at the expected rate

    // Displacement (zero - testing source term only)
    static double disp_x(const double x,
                         const double y,
                         const double t) {
        return 0.0;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double t) {
        return 0.0;
    } // disp_y

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    // Velocity (zero)
    static double vel_x(const double x,
                        const double y,
                        const double t) {
        return 0.0;
    } // vel_x

    static double vel_y(const double x,
                        const double y,
                        const double t) {
        return 0.0;
    } // vel_y

    static const char* vel_units(void) {
        return "m/s";
    } // vel_units

    // Acceleration (zero)
    static double acc_x(const double x,
                        const double y,
                        const double t) {
        return 0.0;
    } // acc_x

    static double acc_y(const double x,
                        const double y,
                        const double t) {
        return 0.0;
    } // acc_y

    static const char* acc_units(void) {
        return "m/s**2";
    } // acc_units

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

        return PETSC_SUCCESS;
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

        return PETSC_SUCCESS;
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

        return PETSC_SUCCESS;
    } // solnkernel_acc

public:

    static
    TestMomentTensorSource_Data* createData(void) {
        TestMomentTensorSource_Data* data = new TestMomentTensorSource_Data();assert(data);

        data->journalName = "MomentTensorSource2D";
        data->isJacobianLinear = true;
        data->tolerance = 1.0e-6;
        data->allowZeroResidual = true; // Zero solution with source gives non-zero residual

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        // Initialize scales first before using them
        pylith::scales::ElasticityScales::setDynamicElasticity(&data->scales, LENGTH_SCALE, VELOCITY_SCALE);
        data->formulation = pylith::problems::Physics::DYNAMIC;

        data->t = TIME_SNAPSHOT;
        data->dt = 0.05;

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

        data->material.setIdentifier("elasticity");
        data->material.setName("material-id=24");
        data->material.setLabelValue(24);

        // Source information
        data->numSourceAuxSubfields = 3;
        static const char* _sourceAuxSubfields[3] = {"moment_tensor", "time_delay", "center_frequency"};
        data->sourceAuxSubfields = _sourceAuxSubfields;
        static const pylith::topology::Field::Discretization _sourceAuxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // moment_tensor
            pylith::topology::Field::Discretization(0, 1), // time_delay
            pylith::topology::Field::Discretization(0, 1), // center_frequency
        };
        data->sourceAuxDiscretizations = const_cast<pylith::topology::Field::Discretization const*>(_sourceAuxDiscretizations);

        data->sourceAuxDB.addValue("moment_tensor_xx", moment_xx, moment_units());
        data->sourceAuxDB.addValue("moment_tensor_yy", moment_yy, moment_units());
        data->sourceAuxDB.addValue("moment_tensor_xy", moment_xy, moment_units());
        data->sourceAuxDB.addValue("moment_tensor_yx", moment_yx, moment_units());
        data->sourceAuxDB.addValue("time_delay", time_delay, time_delay_units());
        data->sourceAuxDB.addValue("center_frequency", center_frequency, center_frequency_units());
        data->sourceAuxDB.setCoordSys(data->cs);

        data->source.setFormulation(data->formulation);
        data->source.setIdentifier("source");
        data->source.setName("moment-tensor-source");
        data->source.setLabelValue(24);
        data->source.setSubfieldName("velocity"); // Apply source to velocity field for dynamic problems

        // Boundary conditions
        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        pylith::bc::DirichletUserFn* bc = NULL;
        data->bcs.resize(2);
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

}; // _MomentTensorSource2D

// ------------------------------------------------------------------------------------------------
pylith::TestMomentTensorSource_Data*
pylith::MomentTensorSource2D::TriP1(void) {
    TestMomentTensorSource_Data* data = pylith::_MomentTensorSource2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";
    data->useAsciiMesh = true;

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP1


// ------------------------------------------------------------------------------------------------
pylith::TestMomentTensorSource_Data*
pylith::MomentTensorSource2D::TriP2(void) {
    TestMomentTensorSource_Data* data = pylith::_MomentTensorSource2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";
    data->useAsciiMesh = true;

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    static const pylith::topology::Field::Discretization _sourceAuxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // moment_tensor
        pylith::topology::Field::Discretization(0, 2), // time_delay
        pylith::topology::Field::Discretization(0, 2), // center_frequency
    };
    data->sourceAuxDiscretizations = const_cast<pylith::topology::Field::Discretization const*>(_sourceAuxDiscretizations);

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP2


// ------------------------------------------------------------------------------------------------
pylith::TestMomentTensorSource_Data*
pylith::MomentTensorSource2D::QuadQ1(void) {
    TestMomentTensorSource_Data* data = pylith::_MomentTensorSource2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";
    data->useAsciiMesh = true;

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ1


// ------------------------------------------------------------------------------------------------
pylith::TestMomentTensorSource_Data*
pylith::MomentTensorSource2D::QuadQ2(void) {
    TestMomentTensorSource_Data* data = pylith::_MomentTensorSource2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";
    data->useAsciiMesh = true;

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    static const pylith::topology::Field::Discretization _sourceAuxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // moment_tensor
        pylith::topology::Field::Discretization(0, 2), // time_delay
        pylith::topology::Field::Discretization(0, 2), // center_frequency
    };
    data->sourceAuxDiscretizations = const_cast<pylith::topology::Field::Discretization const*>(_sourceAuxDiscretizations);

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ2


// End of file
