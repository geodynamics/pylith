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

#include "Gravity2D.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

namespace pylith {
    class _Gravity2D;
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::_Gravity2D {
private:

    static const double LENGTHSCALE;
    static const double TIMESCALE;
    static const double PRESSURESCALE;
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
        return 1.0e+15;
    } // vp

    static const char* vp_units(void) {
        return "m/s";
    } // vp_units

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
        const double accelerationScale = LENGTHSCALE / (TIMESCALE * TIMESCALE);
        const double densityScale = PRESSURESCALE / (velocityScale * velocityScale);
        return density(x,y) / densityScale * GACC / (accelerationScale) * (YMAX/LENGTHSCALE-y);
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

        data->journalName = "Gravity2D";

        data->isJacobianLinear = true;

        data->normalizer.setLengthScale(LENGTHSCALE);
        data->normalizer.setTimeScale(TIMESCALE);
        data->normalizer.setPressureScale(PRESSURESCALE);
        data->normalizer.computeDensityScale();

        delete data->gravityField;data->gravityField = new spatialdata::spatialdb::GravityField();
        data->gravityField->setGravityDir(0.0, -1.0, 0.0);
        data->gravityField->setGravityAcc(GACC);

        // solnDiscretizations set in derived class.

        data->numAuxSubfields = 4;
        static const char* _auxSubfields[4] = {
            "density",
            "gravitational_acceleration",
            "shear_modulus",
            "bulk_modulus",
        };
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // gravitational_acceleration
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        data->auxDB.addValue("density", density, density_units());
        data->auxDB.addValue("vp", vp, vp_units());
        data->auxDB.addValue("vs", vs, vs_units());
        data->auxDB.setCoordSys(data->cs);

        data->material.useBodyForce(false);
        data->rheology.useReferenceState(false);

        data->material.setIdentifier("incompressibleelasticity");
        data->material.setName("material-id=24");
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

}; // TestIsotropicLinearIncompElasticity2D_Gravity
const double pylith::_Gravity2D::LENGTHSCALE = 1.0e+3;
const double pylith::_Gravity2D::TIMESCALE = 2.0;
const double pylith::_Gravity2D::PRESSURESCALE = 2.25e+10;
const double pylith::_Gravity2D::GACC = 9.80665;
const double pylith::_Gravity2D::YMAX = +4.0e+3;

// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::Gravity2D::TriP1(void) {
    TestIncompressibleElasticity_Data* data = pylith::_Gravity2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 1), // density
        pylith::topology::Field::Discretization(0, 1), // gravitational_acceleration
        pylith::topology::Field::Discretization(0, 1), // shear_modulus
        pylith::topology::Field::Discretization(0, 1), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP1


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::Gravity2D::TriP2(void) {
    TestIncompressibleElasticity_Data* data = pylith::_Gravity2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(1, 2), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // gravitational_acceleration
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP2


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::Gravity2D::TriP3(void) {
    TestIncompressibleElasticity_Data* data = pylith::_Gravity2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(2, 3), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // gravitational_acceleration
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP3


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::Gravity2D::QuadQ2(void) {
    TestIncompressibleElasticity_Data* data = pylith::_Gravity2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(1, 2), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // gravitational_acceleration
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ2


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::Gravity2D::QuadQ3(void) {
    TestIncompressibleElasticity_Data* data = pylith::_Gravity2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(2, 3), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // gravitational_acceleration
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ3


// End of file
