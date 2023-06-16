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

#include "QuadTrig.hh" // Implementation of test data

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    class _QuadTrig;
} // pylith

class pylith::_QuadTrig {
    static const double LENGTH_SCALE;
    static const double TIME_SCALE;
    static const double PRESSURE_SCALE;
    static const double TIME_SNAPSHOT; // nondimensional
    static const double AMPLITUDE;

    // Porosity: phi
    static double porosity(const double x,
                          const double y) {
        return 0.1;
    }

    static const char* porosity_units(void) {
        return " ";
    }

    // Solid density
    static double solid_density(const double x,
                          const double y) {
        return 2.5;
    }

    // Fluid density
    static double fluid_density(const double x,
                          const double y) {
        return 1.0;
    }

    // Bulk density: \rho_b = \phi * \rho_f + (1-\phi) * \rho_s
    static double bulk_density(const double x,
                          const double y) {
        const double rho_f = fluid_density(x,y);
        const double rho_s = solid_density(x,y);
        const double phi = porosity(x,y);
        return phi * rho_f + (1-phi) * rho_s;
    }

    // Density units
    static const char* density_units(void) {
        return "kg/m**3";
    }

    // Fluid viscosity: mu_f
    static double fluid_viscosity(const double x,
                     const double y) {
        return 1.0;
    }

    // Fluid viscosity units
    static const char* fluid_viscosity_units(void) {
        return "Pa*s";
    }

    // isotropic_permeability: k
    static double isotropic_permeability(const double x,
                     const double y) {
        return 1.5;
    }

    // isotropic_permeability units
    static const char* isotropic_permeability_units(void) {
        return "m**2";
    }

    // Biot coefficient: alpha
    static double biot_coefficient(const double x,
                     const double y) {
        return 0.6;
    }

    static const char* biot_coefficient_units(void) {
        return " ";
    }


    // Shear Modulus: mu
    static double shear_modulus(const double x,
                     const double y) {
        return 3.0;
    }

    // Drained bulk modulus: K_d
    static double drained_bulk_modulus(const double x,
                     const double y) {
        return 4.0;
    }

    // Fluid bulk modulus: K_f
    static double fluid_bulk_modulus(const double x,
                     const double y) {
        return 8.0;
    }

    // Solid bulk modulus: K_s
    static double solid_bulk_modulus(const double x,
                     const double y) {
        return 10.0;
    }

    // Biot Modulus: M
    static double biot_modulus(const double x,
                     const double y) {
        const double K_f = fluid_bulk_modulus(x,y);
        const double K_s = solid_bulk_modulus(x,y);
        const double phi = porosity(x,y);
        const double alpha = biot_coefficient(x,y);
        const double M_inv = (alpha - phi) / K_s + phi / K_f;
        return 1.0 / M_inv ;
    }

    // Biot Modulus unit
    static const char* modulus_units(void) {
        return "Pa";
    }

    // Solution subfields.

    // Displacement
    static double disp_x(const double x,
                         const double y,
                         const double t) {
        return AMPLITUDE*(x*x)*cos(t);
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double t) {
        return AMPLITUDE*(y*y -2*x*y)*cos(t);
    } // disp_y

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    // Velocity
    static double vel_x(const double x,
                        const double y,
                        const double t) {

        return -AMPLITUDE*(x*x)*sin(t);
    } // vel_x

    static double vel_y(const double x,
                        const double y,
                        const double t) {
        return -AMPLITUDE*(y*y -2*x*y)*sin(t);
    } // vel_y

    static const char* vel_units(void) {
        return "m/s";
    } // vel_units

    // Acceleration
    static double acc_x(const double x,
                        const double y,
                        const double t) {
        return -AMPLITUDE*(x*x)*cos(t);
    } // vel_x

    static double acc_y(const double x,
                        const double y,
                        const double t) {
        return -AMPLITUDE*(y*y -2*x*y)*cos(t);
    } // vel_y

    static const char* acc_units(void) {
        return "m/s**2";
    } // vel_units

    // Pressure
    static double pressure(const double x,
                           const double y,
                           const double t) {
        return (x+y)*cos(t);
    } // pressure

    static double pressure_dot(const double x,
                           const double y,
                           const double t) {
        return -(x+y)*sin(t);
    } // pressure

    static const char* pressure_units(void) {
        return "Pa";
    } // vp_units

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

        s[0] = pressure(x[0], x[1], t);

        return 0;
    } // solnkernel_pressure

    static PetscErrorCode solnkernel_pressure_dot(PetscInt spaceDim,
                                         PetscReal t,
                                         const PetscReal x[],
                                         PetscInt numComponents,
                                         PetscScalar* s,
                                         void* context) {
        assert(2 == spaceDim);
        assert(x);
        assert(1 == numComponents);
        assert(s);

        s[0] = pressure_dot(x[0], x[1], t);

        return 0;
    } // solnkernel_pressure_dot

    // Body force for the momentum equation
    static void fu_body_force(const PylithInt dim,
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
                   PylithScalar fu[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        const double alpha = biot_coefficient(x[0], x[1]);
        const double mu = shear_modulus(x[0], x[1]);
        const double K_d = drained_bulk_modulus(x[0], x[1]);
        const double lambda = K_d - 2.0 * mu / 3.0;
        const double rho_b = bulk_density(x[0], x[1]);

        fu[0] = (-rho_b*x[0]*x[0] + alpha - 2*mu) * cos(t);
        fu[1] = (rho_b*(2*x[0]*x[1] - x[1]*x[1]) + alpha - 2*lambda - 4*mu) * cos(t);
    } // fu_body_force

// Source force for fluid part
    static void fp_source_force(const PylithInt dim,
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
                   PylithScalar fp[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        const double alpha = biot_coefficient(x[0], x[1]);
        const double mu_f = fluid_viscosity(x[0], x[1]);
        const double k = isotropic_permeability(x[0], x[1]);
        const double M = biot_modulus(x[0], x[1]);
        const double rho_f = fluid_density(x[0], x[1]);

        fp[0] = 2.0*k*x[1]*rho_f*cos(t) / mu_f;
        fp[0] -= 2.0*x[1]*alpha*sin(t);
        fp[0] -= (x[0] + x[1])*sin(t) / M;
    } // fp_source_force

public:

    static
    TestLinearPoroElasticity_Data* createData(void) {
        TestLinearPoroElasticity_Data* data = new TestLinearPoroElasticity_Data();assert(data);

        data->journalName = "QuadTrig";
        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        data->normalizer.setLengthScale(LENGTH_SCALE);
        data->normalizer.setTimeScale(TIME_SCALE);
        data->normalizer.setPressureScale(PRESSURE_SCALE);
        data->normalizer.computeDensityScale();
        data->formulation = pylith::problems::Physics::DYNAMIC;

        data->t = TIME_SNAPSHOT;
        data->dt = 0.05;

        // solnDiscretizations set in derived class.

        // Material information
        data->numAuxSubfields = 4;
        static const char* _auxSubfields[4] = {"solid_density", "fluid_density", "fluid_viscosity", "porocity"};
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 1), // solid_density
            pylith::topology::Field::Discretization(0, 1), // fluid_density
            pylith::topology::Field::Discretization(0, 1), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 1), // porocity
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization const*>(_auxDiscretizations);

        data->auxDB.addValue("solid_density", solid_density, density_units());
        data->auxDB.addValue("fluid_density", fluid_density, density_units());
        data->auxDB.addValue("fluid_viscosity", fluid_viscosity, fluid_viscosity_units());
        data->auxDB.addValue("porosity", porosity, porosity_units());
        data->auxDB.addValue("shear_modulus", shear_modulus, modulus_units());
        data->auxDB.addValue("drained_bulk_modulus", drained_bulk_modulus, modulus_units());
        data->auxDB.addValue("biot_coefficient", biot_coefficient, biot_coefficient_units());
        data->auxDB.addValue("fluid_bulk_modulus", fluid_bulk_modulus, modulus_units());
        data->auxDB.addValue("solid_bulk_modulus", solid_bulk_modulus, modulus_units());
        data->auxDB.addValue("isotropic_permeability", isotropic_permeability, isotropic_permeability_units());

        data->auxDB.setCoordSys(data->cs);

        data->material.setFormulation(data->formulation);
        data->material.useBodyForce(false);
        data->rheology.useReferenceState(false);

        data->material.setDescription("Linear PoroElastodynamics");
        data->material.setLabelValue(24);

        // Boundary conditions
        pylith::bc::DirichletUserFn* bc = NULL;
        data->bcs.resize(3);
        {
            static const PylithInt constrainedDOF[2] = {0, 1};
            static const PylithInt numConstrained = 2;
            bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(solnkernel_disp);
            bc->setUserFnDot(solnkernel_vel);
            data->bcs[0] = bc;
        }

        {
            static const PylithInt constrainedDOF[1] = {0};
            static const PylithInt numConstrained = 1;
            bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(solnkernel_pressure);
            bc->setUserFnDot(solnkernel_pressure_dot);
            data->bcs[1] = bc;
        }

        {
            static const PylithInt constrainedDOF[2] = {0, 1};
            static const PylithInt numConstrained = 2;
            bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("velocity");
            bc->setLabelName("boundary");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(solnkernel_vel);
            bc->setUserFnDot(solnkernel_acc);
            data->bcs[2] = bc;
        }

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[3] = {
            solnkernel_disp,
            solnkernel_pressure,
            solnkernel_vel,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        static const pylith::testing::MMSTest::solution_fn _exactSolnDotFns[3] = {
            solnkernel_vel,
            solnkernel_pressure_dot,
            solnkernel_acc,
        };
        data->exactSolnDotFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnDotFns);

        return data;
    } // createData

}; // QuadTrig

const double pylith::_QuadTrig::LENGTH_SCALE = 1.0e+3;
const double pylith::_QuadTrig::TIME_SCALE = 10.0;
const double pylith::_QuadTrig::PRESSURE_SCALE = 3.0e+10;
const double pylith::_QuadTrig::TIME_SNAPSHOT = 7.657345769747113;
const double pylith::_QuadTrig::AMPLITUDE = 1.0e+2;

// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroElasticity_Data*
pylith::QuadTrig::TriP2(void) {
    TestLinearPoroElasticity_Data* data = pylith::_QuadTrig::createData();assert(data);

    data->meshFilename = "data/tri.mesh";
    // data->useAsciiMesh = false;
    // data->tolerance = 2.0e-4;

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(1, 2), // pressure
        pylith::topology::Field::Discretization(2, 2), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP2


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroElasticity_Data*
pylith::QuadTrig::TriP3(void) {
    TestLinearPoroElasticity_Data* data = pylith::_QuadTrig::createData();assert(data);

    data->meshFilename = "data/tri.mesh";
    // data->tolerance = 5.0e-7;

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 3), // solid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 3), // porocity
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(2, 3), // pressure
        pylith::topology::Field::Discretization(3, 3), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP3


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroElasticity_Data*
pylith::QuadTrig::TriP4(void) {
    TestLinearPoroElasticity_Data* data = pylith::_QuadTrig::createData();assert(data);

    data->meshFilename = "data/tri.mesh";
    // data->useAsciiMesh = false;
    // data->tolerance = 1.0e-9;

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 4), // solid_density
        pylith::topology::Field::Discretization(0, 4), // fluid_density
        pylith::topology::Field::Discretization(0, 4), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 4), // porocity
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(4, 4), // disp
        pylith::topology::Field::Discretization(3, 4), // pressure
        pylith::topology::Field::Discretization(4, 4), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP4


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroElasticity_Data*
pylith::QuadTrig::TriP5(void) {
    TestLinearPoroElasticity_Data* data = pylith::_QuadTrig::createData();assert(data);

    data->meshFilename = "data/tri.mesh";
    // data->tolerance = 1.0e-9;

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 5), // solid_density
        pylith::topology::Field::Discretization(0, 5), // fluid_density
        pylith::topology::Field::Discretization(0, 5), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 5), // porocity
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(5, 5), // disp
        pylith::topology::Field::Discretization(4, 5), // pressure
        pylith::topology::Field::Discretization(5, 5), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP5


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroElasticity_Data*
pylith::QuadTrig::QuadQ2(void) {
    TestLinearPoroElasticity_Data* data = pylith::_QuadTrig::createData();assert(data);

    data->meshFilename = "data/quad.mesh";
    // data->tolerance = 1.0e-4;

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(1, 2), // pressure
        pylith::topology::Field::Discretization(2, 2), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ2


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroElasticity_Data*
pylith::QuadTrig::QuadQ2Distorted(void) {
    TestLinearPoroElasticity_Data* data = pylith::_QuadTrig::createData();assert(data);

    data->meshFilename = "data/quad_distorted.mesh";
    // data->tolerance = 1.0e-5;

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(1, 2), // pressure
        pylith::topology::Field::Discretization(2, 2), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ2Distorted


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroElasticity_Data*
pylith::QuadTrig::QuadQ3(void) {
    TestLinearPoroElasticity_Data* data = pylith::_QuadTrig::createData();assert(data);

    data->meshFilename = "data/quad.mesh";
    // data->useAsciiMesh = false;
    // data->tolerance = 2.0e-8;

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 3), // solid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 3), // porocity
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(2, 3), // pressure
        pylith::topology::Field::Discretization(3, 3), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ3


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroElasticity_Data*
pylith::QuadTrig::QuadQ4(void) {
    TestLinearPoroElasticity_Data* data = pylith::_QuadTrig::createData();assert(data);

    data->meshFilename = "data/quad.mesh";
    // data->tolerance = 1.0e-9;

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 4), // solid_density
        pylith::topology::Field::Discretization(0, 4), // fluid_density
        pylith::topology::Field::Discretization(0, 4), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 4), // porocity
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(4, 4), // disp
        pylith::topology::Field::Discretization(3, 4), // pressure
        pylith::topology::Field::Discretization(4, 4), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ4


// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroElasticity_Data*
pylith::QuadTrig::QuadQ5(void) {
    TestLinearPoroElasticity_Data* data = pylith::_QuadTrig::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 5), // solid_density
        pylith::topology::Field::Discretization(0, 5), // fluid_density
        pylith::topology::Field::Discretization(0, 5), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 5), // porocity
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(5, 5), // disp
        pylith::topology::Field::Discretization(4, 5), // disp
        pylith::topology::Field::Discretization(5, 5), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ5


// End of file
