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

#include "PressureGradient.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

namespace pylith {
    class _PressureGradient;
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::_PressureGradient {
    static const double LENGTH_SCALE;
    static const double TIME_SCALE;
    static const double PRESSURE_SCALE;

    static const double PRESSURE0; // dimensional
    static const double XMAX; // dimensional

    // Density
    static double solid_density(const double x,
                                const double y) {
        return 2500.0;
    } // solid_density

    static double fluid_density(const double x,
                                const double y) {
        return 1000.0;
    } // fluid_density

    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    // Fluid viscosity
    static double fluid_viscosity(const double x,
                                  const double y) {
        return 1.0e-3;
    } // vs

    static const char* viscosity_units(void) {
        return "Pa*s";
    } // fluid_viscosity_units

    // Porosity
    static double porosity(const double x,
                           const double y) {
        return 0.02;
    } // porosity

    static const char* porosity_units(void) {
        return "none";
    } // porosity_units

    // Shear modulus
    static double shear_modulus(const double x,
                                const double y) {
        return 3.0e+10;
    } // shear_modulus

    static const char* modulus_units(void) {
        return "Pa";
    } // modulus_units

    // Drained bulk modulus
    static double drained_bulk_modulus(const double x,
                                       const double y) {
        return 8.0e+10;
    } // drained_bulk_modulus

    // Biot coefficient
    static double biot_coefficient(const double x,
                                   const double y) {
        return 0.8;
    } // biot_coefficient

    static const char* biot_coefficient_units(void) {
        return "none";
    } // biot_coefficient_units

    // Fluid modulus
    static double fluid_bulk_modulus(const double x,
                                     const double y) {
        return 1.0e+10;
    } // fluid_bulk_modulus

    // Permeability
    static double isotropic_permeability(const double x,
                                         const double y) {
        return 1.0e-14;
    } // isotropic_permeability

    static const char* permeability_units(void) {
        return "m*m";
    } // permeability_units

    // Solution subfields (nondimensional)

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        const double muN = shear_modulus(x, y) / PRESSURE_SCALE;
        const double lambdaN = drained_bulk_modulus(x, y) / PRESSURE_SCALE - 2.0/3.0 * muN;
        const double alpha = biot_coefficient(x, y);
        return -0.5 * alpha  * (PRESSURE0 / PRESSURE_SCALE) / (lambdaN + 2.0*muN) * (x*x / (XMAX / LENGTH_SCALE));
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return 0.0;
    } // disp_y

    // Pressure
    static double fluid_pressure(const double x,
                                 const double y) {
        return (PRESSURE0 / PRESSURE_SCALE) * (1.0 - x / (XMAX / LENGTH_SCALE));
    } // fluid_pressure

    // Trace strain
    static double trace_strain(const double x,
                               const double y) {
        const double muN = shear_modulus(x, y) / PRESSURE_SCALE;
        const double lambdaN = drained_bulk_modulus(x, y) / PRESSURE_SCALE - 2.0/3.0 * muN;
        const double alpha = biot_coefficient(x, y);
        return -alpha  * (PRESSURE0 / PRESSURE_SCALE)  / (lambdaN + 2.0*muN) * (x / (XMAX / LENGTH_SCALE));
    } // trace_strain

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

    static PetscErrorCode solnkernel_fluid_pressure(PetscInt spaceDim,
                                                    PetscReal t,
                                                    const PetscReal x[],
                                                    PetscInt numComponents,
                                                    PetscScalar* s,
                                                    void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = fluid_pressure(x[0], x[1]);

        return 0;
    } // solnkernel_fluid_pressure

    static PetscErrorCode solnkernel_trace_strain(PetscInt spaceDim,
                                                  PetscReal t,
                                                  const PetscReal x[],
                                                  PetscInt numComponents,
                                                  PetscScalar* s,
                                                  void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = trace_strain(x[0], x[1]);

        return 0;
    } // solnkernel_trace_strain

    static PetscErrorCode solnkernel_vel(PetscInt spaceDim,
                                         PetscReal t,
                                         const PetscReal x[],
                                         PetscInt numComponents,
                                         PetscScalar* s,
                                         void* context) {
        assert(2 == spaceDim);
        assert(2 == numComponents);
        assert(s);

        s[0] = 0.0;
        s[1] = 0.0;

        return 0;
    } // solnkernel_vel

    static PetscErrorCode solnkernel_fluid_pressure_dot(PetscInt spaceDim,
                                                        PetscReal t,
                                                        const PetscReal x[],
                                                        PetscInt numComponents,
                                                        PetscScalar* s,
                                                        void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = 0.0;

        return 0;
    } // solnkernel_fluid_pressure_dot

    static PetscErrorCode solnkernel_trace_strain_dot(PetscInt spaceDim,
                                                      PetscReal t,
                                                      const PetscReal x[],
                                                      PetscInt numComponents,
                                                      PetscScalar* s,
                                                      void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = 0.0;

        return 0;
    } // solnkernel_trace_strain_dot

public:

    static
    TestLinearPoroelasticity_Data* createData(void) {
        TestLinearPoroelasticity_Data* data = new TestLinearPoroelasticity_Data();assert(data);

        data->journalName = "PressureGradient";
        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        data->normalizer.setLengthScale(LENGTH_SCALE);
        data->normalizer.setTimeScale(TIME_SCALE);
        data->normalizer.setPressureScale(PRESSURE_SCALE);
        data->normalizer.computeDensityScale();

        // solnDiscretizations set in derived class.

        // Material information
        data->numAuxSubfields = 9;
        static const char* _auxSubfields[9] = { // order must match order of subfields in auxiliary field
            "solid_density",
            "fluid_density",
            "fluid_viscosity",
            "porosity",
            "shear_modulus",
            "drained_bulk_modulus",
            "biot_coefficient",
            "biot_modulus",
            "isotropic_permeability",
        };
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
            pylith::topology::Field::Discretization(0, 1), // solid_density
            pylith::topology::Field::Discretization(0, 1), // fluid_density
            pylith::topology::Field::Discretization(0, 1), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 1), // porosity
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // drained_bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // biot_coefficient
            pylith::topology::Field::Discretization(0, 1), // biot_modulus
            pylith::topology::Field::Discretization(0, 1), // isotropic_permeability
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        data->auxDB.addValue("solid_density", solid_density, density_units());
        data->auxDB.addValue("fluid_density", fluid_density, density_units());
        data->auxDB.addValue("fluid_viscosity", fluid_viscosity, viscosity_units());
        data->auxDB.addValue("porosity", porosity, porosity_units());
        data->auxDB.addValue("shear_modulus", shear_modulus, modulus_units());
        data->auxDB.addValue("drained_bulk_modulus", drained_bulk_modulus, modulus_units());
        data->auxDB.addValue("biot_coefficient", biot_coefficient, modulus_units());
        data->auxDB.addValue("fluid_bulk_modulus", fluid_bulk_modulus, modulus_units());
        data->auxDB.addValue("isotropic_permeability", isotropic_permeability, permeability_units());
        data->auxDB.setCoordSys(data->cs);

        data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
        data->rheology.useReferenceState(false);

        data->material.setIdentifier("poroelasticity");
        data->material.setName("material-id=24");
        data->material.setLabelValue(24);

        static const PylithInt constrainedX[1] = { 0 };
        static const PylithInt constrainedY[1] = { 1 };
        static const PylithInt numConstrained = 1;
        data->bcs.resize(6);
        { // Displacement -x
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedX, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[0] = bc;
        }
        { // Displacement +x
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedX, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[1] = bc;
        }
        { // Displacement -y
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_yneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedY, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[2] = bc;
        }
        { // Displacement +y
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_ypos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedY, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[3] = bc;
        }
        { // Pressure -x
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedX, numConstrained);
            bc->setUserFn(solnkernel_fluid_pressure);
            data->bcs[4] = bc;
        }
        { // Pressure +x
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedX, numConstrained);
            bc->setUserFn(solnkernel_fluid_pressure);
            data->bcs[5] = bc;
        }

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[3] = {
            solnkernel_disp,
            solnkernel_fluid_pressure,
            solnkernel_trace_strain,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        data->exactSolnDotFns = nullptr;

        return data;
    } // createData

    static
    TestLinearPoroelasticity_Data* createDataStateVars(void) {
        TestLinearPoroelasticity_Data* data = createData();

        data->material.useStateVars(true);

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[6] = {
            solnkernel_disp,
            solnkernel_fluid_pressure,
            solnkernel_trace_strain,
            solnkernel_vel,
            solnkernel_fluid_pressure_dot,
            solnkernel_trace_strain_dot,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);

        return data;
    } // createDataStateVars

}; // PressureGradient
const double pylith::_PressureGradient::LENGTH_SCALE = 1.0e+3;
const double pylith::_PressureGradient::TIME_SCALE = 2.0;
const double pylith::_PressureGradient::PRESSURE_SCALE = 2.25e+10;

const double pylith::_PressureGradient::PRESSURE0 = 4.0e+6;
const double pylith::_PressureGradient::XMAX = 8.0e+3;

// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroelasticity_Data*
pylith::PressureGradient::TriP2P1P1(void) {
    TestLinearPoroelasticity_Data* data = pylith::_PressureGradient::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(1, 2), // fluid pressure
        pylith::topology::Field::Discretization(1, 2), // trace strain
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
        pylith::topology::Field::Discretization(0, 2), // solid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 2), // porosity
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 2), // biot_coefficient
        pylith::topology::Field::Discretization(0, 2), // biot_modulus
        pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP2P1P1


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroelasticity_Data*
pylith::PressureGradient::TriP3P2P2(void) {
    TestLinearPoroelasticity_Data* data = pylith::_PressureGradient::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(3, 3), // displacement
        pylith::topology::Field::Discretization(2, 3), // fluid pressure
        pylith::topology::Field::Discretization(2, 3), // trace strain
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
        pylith::topology::Field::Discretization(0, 3), // solid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 3), // porosity
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 3), // biot_coefficient
        pylith::topology::Field::Discretization(0, 3), // biot_modulus
        pylith::topology::Field::Discretization(0, 3), // isotropic_permeability
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP3P2P2


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroelasticity_Data*
pylith::PressureGradient::QuadQ2Q1Q1(void) {
    TestLinearPoroelasticity_Data* data = pylith::_PressureGradient::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(1, 2), // fluid pressure
        pylith::topology::Field::Discretization(1, 2), // trace strain
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
        pylith::topology::Field::Discretization(0, 2), // solid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 2), // porosity
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 2), // biot_coefficient
        pylith::topology::Field::Discretization(0, 2), // biot_modulus
        pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ2Q1Q1


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroelasticity_Data*
pylith::PressureGradient::QuadQ3Q2Q2(void) {
    TestLinearPoroelasticity_Data* data = pylith::_PressureGradient::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(3, 3), // displacement
        pylith::topology::Field::Discretization(2, 3), // fluid pressure
        pylith::topology::Field::Discretization(2, 3), // trace strain
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
        pylith::topology::Field::Discretization(0, 3), // solid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 3), // porosity
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 3), // biot_coefficient
        pylith::topology::Field::Discretization(0, 3), // biot_modulus
        pylith::topology::Field::Discretization(0, 3), // isotropic_permeability
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ2Q1Q1


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroelasticity_Data*
pylith::PressureGradient::TriP2P1P1_StateVars(void) {
    TestLinearPoroelasticity_Data* data = pylith::_PressureGradient::createDataStateVars();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 6;
    static const pylith::topology::Field::Discretization _solnDiscretizations[6] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(1, 2), // fluid pressure
        pylith::topology::Field::Discretization(1, 2), // trace strain
        pylith::topology::Field::Discretization(2, 2), // velocity
        pylith::topology::Field::Discretization(1, 2), // fluid pressure dot
        pylith::topology::Field::Discretization(1, 2), // trace strain dot
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
        pylith::topology::Field::Discretization(0, 2), // solid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 2), // porosity
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 2), // biot_coefficient
        pylith::topology::Field::Discretization(0, 2), // biot_modulus
        pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP2P1P1_StateVars


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroelasticity_Data*
pylith::PressureGradient::TriP3P2P2_StateVars(void) {
    TestLinearPoroelasticity_Data* data = pylith::_PressureGradient::createDataStateVars();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 6;
    static const pylith::topology::Field::Discretization _solnDiscretizations[6] = {
        pylith::topology::Field::Discretization(3, 3), // displacement
        pylith::topology::Field::Discretization(2, 3), // fluid pressure
        pylith::topology::Field::Discretization(2, 3), // trace strain
        pylith::topology::Field::Discretization(3, 3), // velocity
        pylith::topology::Field::Discretization(2, 3), // fluid pressure dot
        pylith::topology::Field::Discretization(2, 3), // trace strain dot
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
        pylith::topology::Field::Discretization(0, 3), // solid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 3), // porosity
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 3), // biot_coefficient
        pylith::topology::Field::Discretization(0, 3), // biot_modulus
        pylith::topology::Field::Discretization(0, 3), // isotropic_permeability
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP3P2P2_StateVars


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroelasticity_Data*
pylith::PressureGradient::QuadQ2Q1Q1_StateVars(void) {
    TestLinearPoroelasticity_Data* data = pylith::_PressureGradient::createDataStateVars();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 6;
    static const pylith::topology::Field::Discretization _solnDiscretizations[6] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(1, 2), // fluid pressure
        pylith::topology::Field::Discretization(1, 2), // trace strain
        pylith::topology::Field::Discretization(2, 2), // velocity
        pylith::topology::Field::Discretization(1, 2), // fluid pressure dot
        pylith::topology::Field::Discretization(1, 2), // trace strain dot
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
        pylith::topology::Field::Discretization(0, 2), // solid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 2), // porosity
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 2), // biot_coefficient
        pylith::topology::Field::Discretization(0, 2), // biot_modulus
        pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ2Q1Q1_StateVars


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroelasticity_Data*
pylith::PressureGradient::QuadQ3Q2Q2_StateVars(void) {
    TestLinearPoroelasticity_Data* data = pylith::_PressureGradient::createDataStateVars();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 6;
    static const pylith::topology::Field::Discretization _solnDiscretizations[6] = {
        pylith::topology::Field::Discretization(3, 3), // displacement
        pylith::topology::Field::Discretization(2, 3), // fluid pressure
        pylith::topology::Field::Discretization(2, 3), // trace strain
        pylith::topology::Field::Discretization(3, 3), // velocity
        pylith::topology::Field::Discretization(2, 3), // fluid pressure dot
        pylith::topology::Field::Discretization(2, 3), // trace strain dot
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
        pylith::topology::Field::Discretization(0, 3), // solid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 3), // porosity
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 3), // biot_coefficient
        pylith::topology::Field::Discretization(0, 3), // biot_modulus
        pylith::topology::Field::Discretization(0, 3), // isotropic_permeability
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ3Q2Q2_StateVars


// End of file
