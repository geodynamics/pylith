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

#include "BodyForce2D.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

namespace pylith {
    class _BodyForce2D;
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::_BodyForce2D {
private:

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
        return 1.0e+15;
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

    // Solution fields (nondimensional)

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return 0.0 / LENGTHSCALE;
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return 0.0 / LENGTHSCALE;
    } // disp_y

    // Pressure
    static double pressure(const double x,
                           const double y) {
        const double velocityScale = LENGTHSCALE / TIMESCALE;
        const double densityScale = PRESSURESCALE / (velocityScale * velocityScale);
        const double accelerationScale = LENGTHSCALE / (TIMESCALE * TIMESCALE);
        const double forceScale = densityScale * accelerationScale;
        const double bodyforceN = BODYFORCE / forceScale;
        return -bodyforceN * (XMAX/LENGTHSCALE - x);
    } // pressure

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        assert(2 == spaceDim);
        assert(2 == numComponents);
        assert(s);
        assert(x);

        s[0] = disp_x(x[0], x[1]);
        s[1] = disp_y(x[0], x[1]);

        return PETSC_SUCCESS;
    } // solnkernel_disp

    static PetscErrorCode solnkernel_pressure(PetscInt spaceDim,
                                              PetscReal t,
                                              const PetscReal x[],
                                              PetscInt numComponents,
                                              PetscScalar* s,
                                              void* context) {
        assert(2 == spaceDim);
        assert(x);
        assert(1 == numComponents);
        assert(s);

        s[0] = pressure(x[0], x[1]);

        return PETSC_SUCCESS;
    } // solnkernel_pressure

public:

    static
    TestIncompressibleElasticity_Data* createData(void) {
        TestIncompressibleElasticity_Data* data = new TestIncompressibleElasticity_Data();assert(data);

        data->journalName = "BodyForce2D";

        data->isJacobianLinear = true;

        data->normalizer.setLengthScale(LENGTHSCALE);
        data->normalizer.setTimeScale(TIMESCALE);
        data->normalizer.setPressureScale(PRESSURESCALE);
        data->normalizer.computeDensityScale();

        // solnDiscretizations set in derived class.

        data->numAuxSubfields = 4;
        static const char* _auxSubfields[4] = {
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

        data->material.useBodyForce(true);
        data->rheology.useReferenceState(false);

        data->material.setDescription("Isotropic Linear Incompressible Elasticity Plane Strain");
        data->material.setLabelValue(24);

        data->bcs.resize(2);
        { // disp
            static const PylithInt constrainedDOF[2] = {0, 1};
            static const PylithInt numConstrainedDOF = 2;
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrainedDOF);
            bc->setUserFn(solnkernel_disp);
            data->bcs[0] = bc;
        } // disp

        { // pressure
            static const PylithInt constrainedDOF[1] = {0};
            static const PylithInt numConstrainedDOF = 1;
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrainedDOF);
            bc->setUserFn(solnkernel_pressure);
            data->bcs[1] = bc;
        } // pressure

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[2] = {
            solnkernel_disp,
            solnkernel_pressure,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        data->exactSolnDotFns = nullptr;

        return data;
    } // createData

}; // TestIsotropicLinearIncompElasticity2D_BodyForce
const double pylith::_BodyForce2D::LENGTHSCALE = 1.0e+3;
const double pylith::_BodyForce2D::TIMESCALE = 2.0;
const double pylith::_BodyForce2D::PRESSURESCALE = 2.25e+10;
const double pylith::_BodyForce2D::BODYFORCE = 20.0e+3;
const double pylith::_BodyForce2D::XMAX = +4.0e+3;

// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::BodyForce2D::TriP1(void) {
    TestIncompressibleElasticity_Data* data = pylith::_BodyForce2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 1), // density
        pylith::topology::Field::Discretization(0, 1), // body_force
        pylith::topology::Field::Discretization(0, 1), // shear_modulus
        pylith::topology::Field::Discretization(0, 1), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP1


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::BodyForce2D::TriP2(void) {
    TestIncompressibleElasticity_Data* data = pylith::_BodyForce2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(1, 2), // pressure
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
pylith::TestIncompressibleElasticity_Data*
pylith::BodyForce2D::TriP3(void) {
    TestIncompressibleElasticity_Data* data = pylith::_BodyForce2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(2, 3), // pressure
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
pylith::TestIncompressibleElasticity_Data*
pylith::BodyForce2D::QuadQ2(void) {
    TestIncompressibleElasticity_Data* data = pylith::_BodyForce2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(1, 2), // pressure
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
} // QuadQ2


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::BodyForce2D::QuadQ3(void) {
    TestIncompressibleElasticity_Data* data = pylith::_BodyForce2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(2, 3), // pressure
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
