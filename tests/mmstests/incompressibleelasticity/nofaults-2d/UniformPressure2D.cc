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

#include "UniformPressure2D.hh" // ISA UniformPressure2D

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "pylith/scales/ElasticityScales.hh" // USES ElasticityScales

namespace pylith {
    class _UniformPressure2D;
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::_UniformPressure2D {
private:

    static pylith::scales::Scales scales;
    static const double PRESSURE; // nondimensional

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
        return 1.0e+20;
    } // vp

    static const char* vp_units(void) {
        return "m/s";
    } // vp_units

    // Solution subfields.

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        const double rigidityScale = scales.getRigidityScale();
        const double bulkModulus = density(x, y) * (vp(x,y)*vp(x,y) - 4.0/3.0*vs(x,y)*vs(x,y)) / rigidityScale;

        return -0.5 * PRESSURE / bulkModulus * x;
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        const double rigidityScale = scales.getRigidityScale();
        const double bulkModulus = density(x, y) * (vp(x,y)*vp(x,y) - 4.0/3.0*vs(x,y)*vs(x,y)) / rigidityScale;

        return -0.5 * PRESSURE / bulkModulus * y;
    } // disp_y

    // Pressure
    static double pressure(const double x,
                           const double y) {
        return PRESSURE;
    } // pressure

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

        data->journalName = "UniformPressure2D";

        data->isJacobianLinear = true;
        data->jacobianConvergenceRate = 1.0;
        data->tolerance = 2.0e-8;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        scales = data->scales;

        // solnDiscretizations set in derived class.

        data->numAuxSubfields = 3;
        static const char* _auxSubfields[3] = {"density", "shear_modulus", "bulk_modulus"};
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        data->auxDB.addValue("density", density, density_units());
        data->auxDB.addValue("vp", vp, vp_units());
        data->auxDB.addValue("vs", vs, vs_units());
        data->auxDB.setCoordSys(data->cs);

        data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
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

}; // _UniformPressure2D
pylith::scales::Scales pylith::_UniformPressure2D::scales;
const double pylith::_UniformPressure2D::PRESSURE = 3.0e+6;

// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::UniformPressure2D::TriP1(void) {
    TestIncompressibleElasticity_Data* data = pylith::_UniformPressure2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(0, 1), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP1


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::UniformPressure2D::TriP2(void) {
    TestIncompressibleElasticity_Data* data = pylith::_UniformPressure2D::createData();assert(data);

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
        pylith::topology::Field::Discretization(0, 2), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP2


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::UniformPressure2D::TriP3(void) {
    TestIncompressibleElasticity_Data* data = pylith::_UniformPressure2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(0, 3), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP3


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::UniformPressure2D::TriP4(void) {
    TestIncompressibleElasticity_Data* data = pylith::_UniformPressure2D::createData();assert(data);

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
        pylith::topology::Field::Discretization(0, 4), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP4


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::UniformPressure2D::QuadQ1(void) {
    TestIncompressibleElasticity_Data* data = pylith::_UniformPressure2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(0, 1), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ1


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::UniformPressure2D::QuadQ2(void) {
    TestIncompressibleElasticity_Data* data = pylith::_UniformPressure2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(0, 2), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ2


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::UniformPressure2D::QuadQ3(void) {
    TestIncompressibleElasticity_Data* data = pylith::_UniformPressure2D::createData();assert(data);

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
        pylith::topology::Field::Discretization(0, 3), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ3


// ------------------------------------------------------------------------------------------------
pylith::TestIncompressibleElasticity_Data*
pylith::UniformPressure2D::QuadQ4(void) {
    TestIncompressibleElasticity_Data* data = pylith::_UniformPressure2D::createData();assert(data);

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
        pylith::topology::Field::Discretization(0, 4), // pressure
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ4


// End of file
