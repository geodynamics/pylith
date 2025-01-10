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

#include "BodyForce2D.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

namespace pylith {
    class _BodyForce2D;
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::_BodyForce2D {
    static const double LENGTHSCALE;
    static const double TIMESCALE;
    static const double PRESSURESCALE;
    static const double BODYFORCE;
    static const double XMAX;

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

    static double bodyforce_x(const double x,
                              const double y) {
        return BODYFORCE;
    } // bodyforce_x

    static double bodyforce_y(const double x,
                              const double y) {
        return 0.0;
    } // bodyforce_y

    static const char* bodyforce_units(void) {
        return "kg/(m**2*s**2)";
    } // bodyforce_units

    // Solution subfields (nondimensional)

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        const double velocityScale = LENGTHSCALE / TIMESCALE;
        const double densityScale = PRESSURESCALE / (velocityScale * velocityScale);
        const double accelerationScale = LENGTHSCALE / (TIMESCALE * TIMESCALE);
        const double forceScale = densityScale * accelerationScale;
        const double bodyforceN = BODYFORCE / forceScale;
        const double muN = density(x,y) * vs(x,y) * vs(x,y) / PRESSURESCALE;
        const double lambdaN = density(x,y) * vp(x,y) * vp(x,y) / PRESSURESCALE - 2.0*muN;
        const double xp = x - XMAX / LENGTHSCALE;
        return -0.5 * bodyforceN / (lambdaN + 2.0*muN) * (xp*xp);
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
        assert(2 == spaceDim);
        assert(2 == numComponents);
        assert(s);

        s[0] = disp_x(x[0], x[1]);
        s[1] = disp_y(x[0], x[1]);

        return 0;
    } // solnkernel_disp

public:

    static
    TestLinearElasticity_Data* createData(void) {
        TestLinearElasticity_Data* data = new TestLinearElasticity_Data();assert(data);

        data->journalName = "BodyForce2D";
        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        data->normalizer.setLengthScale(LENGTHSCALE);
        data->normalizer.setTimeScale(TIMESCALE);
        data->normalizer.setPressureScale(PRESSURESCALE);
        data->normalizer.computeDensityScale();

        // solnDiscretizations set in derived class.

        // Material information
        data->numAuxSubfields = 4;
        static const char* _auxSubfields[4] = { // order must match order of subfields in auxiliary field
            "density",
            "body_force",
            "shear_modulus",
            "bulk_modulus",
        };
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // body_force
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        data->auxDB.addValue("density", density, density_units());
        data->auxDB.addValue("vp", vp, vp_units());
        data->auxDB.addValue("vs", vs, vs_units());
        data->auxDB.addValue("body_force_x", bodyforce_x, bodyforce_units());
        data->auxDB.addValue("body_force_y", bodyforce_y, bodyforce_units());
        data->auxDB.setCoordSys(data->cs);

        data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
        data->material.useBodyForce(true);
        data->rheology.useReferenceState(false);

        data->material.setDescription("Isotropic Linear Elasticity Plane Strain");
        data->material.setLabelValue(24);

        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        data->bcs.resize(1);
        pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
        bc->setSubfieldName("displacement");
        bc->setLabelName("boundary");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(constrainedDOF, numConstrained);
        bc->setUserFn(solnkernel_disp);
        data->bcs[0] = bc;

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[1] = {
            solnkernel_disp,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        data->exactSolnDotFns = nullptr;

        return data;
    } // createData

}; // BodyForce2D
const double pylith::_BodyForce2D::LENGTHSCALE = 1.0e+3;
const double pylith::_BodyForce2D::TIMESCALE = 2.0;
const double pylith::_BodyForce2D::PRESSURESCALE = 2.25e+10;
const double pylith::_BodyForce2D::BODYFORCE = 5.0e+3;
const double pylith::_BodyForce2D::XMAX = 4.0e+3;

// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::BodyForce2D::TriP2(void) {
    TestLinearElasticity_Data* data = pylith::_BodyForce2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(2, 2), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // body_force
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP2


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::BodyForce2D::TriP3(void) {
    TestLinearElasticity_Data* data = pylith::_BodyForce2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(3, 3), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // body_force
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP3


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::BodyForce2D::QuadQ2(void) {
    TestLinearElasticity_Data* data = pylith::_BodyForce2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(2, 2), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // body_force
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ1


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::BodyForce2D::QuadQ3(void) {
    TestLinearElasticity_Data* data = pylith::_BodyForce2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(3, 3), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // body_force
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ3


// End of file
