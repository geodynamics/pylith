// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestIsotropicLinearPoroelasticity.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh"                  // USES TimeDependent
#include "pylith/materials/Poroelasticity.hh"                // USES Poroelasticity
#include "pylith/materials/IsotropicLinearPoroelasticity.hh" // USES IsotropicLinearPoroelasticity
#include "pylith/fekernels/IsotropicLinearPoroelasticity.hh" // USES IsotropicLinearIncompElasticity kernels
#include "pylith/bc/DirichletUserFn.hh"                      // USES DirichletUserFn

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh"         // USES CSCart
#include "spatialdata/units/Nondimensional.hh"     // USES Nondimensional

namespace pylith
{
    namespace mmstests
    {
        class TestIsotropicLinearPoroelasticity2D_TS_LT;

        class TestIsotropicLinearPoroelasticity2D_TS_LT_TriP2;
        class TestIsotropicLinearPoroelasticity2D_TS_LT_TriP3;
        class TestIsotropicLinearPoroelasticity2D_TS_LT_TriP4;

        class TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ2;
        class TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ3;
        class TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ4;

    } // materials
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT : public pylith::mmstests::TestIsotropicLinearPoroelasticity
{
    static const double LENGTHSCALE;
    static const double TIMESCALE;
    static const double PRESSURESCALE;
    static const double BODYFORCE;
    static const double XMIN;
    static const double XMAX;

    /// Spatial database user functions for auxiiliary subfields (includes derived fields).

    // Porosity

    static double porosity(const double x,
                           const double y)
    {
        return 0.10;
    } // porosity

    static const char *porosity_units(void)
    {
        return "none";
    } // porosity_units

    // Solid Density
    static double solid_density(const double x,
                                const double y)
    {
        return 2500.0;
    } // solid_density

    static const char *solid_density_units(void)
    {
        return "kg/m**3";
    } // solid_density_units

    // Fluid Density
    static double fluid_density(const double x,
                                const double y)
    {
        return 1000.0;
    } // fluid_density

    static const char *fluid_density_units(void)
    {
        return "kg/m**3";
    } // fluid_density_units

    // Fluid viscosity
    static double fluid_viscosity(const double x,
                                  const double y)
    {
        return 1.0;
    } // fluid_viscosity

    static const char *fluid_viscosity_units(void)
    {
        return "Pa*s";
    } // fluid_viscosity_units

    // Shear modulus
    static double shear_modulus(const double x,
                                const double y)
    {
        return 3.0;
    } // shear_modulus

    static const char *shear_modulus_units(void)
    {
        return "Pa";
    } // shear_modulus_units

    // Drained Bulk Modulus
    static double drained_bulk_modulus(const double x,
                                       const double y)
    {
        return 4.0;
    } // drained_bulk_modulus

    static const char *drained_bulk_modulus_units(void)
    {
        return "Pa";
    } // drained_bulk_modulus_units

    // Biot coefficient
    static double biot_coefficient(const double x,
                                   const double y)
    {
        return 0.6;
    } // biot_coefficient

    static const char *biot_coefficient_units(void)
    {
        return "none";
    } // biot_coefficient_units

    // Fluid Bulk Modulus
    static double fluid_bulk_modulus(const double x,
                                     const double y)
    {
        return 8.0;
    } // fluid_bulk_modulus

    static const char *fluid_bulk_modulus_units(void)
    {
        return "Pa";
    } // fluid_bulk_modulus_units

    static double solid_bulk_modulus(const double x,
                                     const double y)
    {
        return 10.0;
    } // solid_bulk_modulus

    static const char *solid_bulk_modulus_units(void)
    {
        return "Pa";
    } // solid_bulk_modulus_units

    // Isotropic permeability
    static double isotropic_permeability(const double x,
                                         const double y)
    {
        return 1.5;
    } // isotropic_permeability

    static const char *isotropic_permeability_units(void)
    {
        return "m**2";
    } // isotropic_permeability_units

    // Derived Fields

    static double biot_modulus(const double x,
                               const double y)
    {
        return 1.0 / (porosity(x, y) / fluid_bulk_modulus(x, y) +
                      (biot_coefficient(x, y) - porosity(x, y)) / solid_bulk_modulus(x, y));
    }

    static const char *biot_modulus_units(void)
    {
        return "Pa";
    } // biot_modulus_units

    // Solution subfields (nondimensional)

    static double displacement_x(const double x,
                                 const double y,
                                 const double t)
    {
        return PetscSinReal(2.0 * PETSC_PI * x);
    } // displacement_x

    static double displacement_y(const double x,
                                 const double y,
                                 const double t)
    {
        return PetscSinReal(2.0 * PETSC_PI * y) - 2.0 * x * y;
    } // displacement_y

    static double pressure(const double x,
                           const double y,
                           const double t)
    {
        return (PetscCosReal(2.0 * PETSC_PI * x) + PetscCosReal(2.0 * PETSC_PI * y)) * t;
    } // pressure

    static double trace_strain(const double x,
                               const double y,
                               const double t)
    {
        return 2.0 * PETSC_PI * (PetscCosReal(2.0 * PETSC_PI * x) + PetscCosReal(2.0 * PETSC_PI * y)) - 2.0 * x;
    } // trace_strain

    static double displacement_t_x(const double x,
                                   const double y,
                                   const double t)
    {
        return 0.0;
    } // displacement_t_x

    static double displacement_t_y(const double x,
                                   const double y,
                                   const double t)
    {
        return 0.0;
    } // displacement_t_y

    static double pressure_t(const double x,
                             const double y,
                             const double t)
    {
        return PetscCosReal(2.0 * PETSC_PI * x) + PetscCosReal(2.0 * PETSC_PI * y);
    } // pressure_t

    static double trace_strain_t(const double x,
                                 const double y,
                                 const double t)
    {
        return 0.0;
    } // trace_strain_t

    static PetscErrorCode solnkernel_displacement(PetscInt spaceDim,
                                                  PetscReal t,
                                                  const PetscReal x[],
                                                  PetscInt numComponents,
                                                  PetscScalar *s,
                                                  void *context)
    {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = displacement_x(x[0], x[1], t);
        s[1] = displacement_y(x[0], x[1], t);

        return 0;
    } // solnkernel_displacement

    static PetscErrorCode solnkernel_pressure(PetscInt spaceDim,
                                              PetscReal t,
                                              const PetscReal x[],
                                              PetscInt numComponents,
                                              PetscScalar *s,
                                              void *context)
    {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(1 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = pressure(x[0], x[1], t);

        return 0;
    } // solnkernel_pressure

    static PetscErrorCode solnkernel_trace_strain(PetscInt spaceDim,
                                                  PetscReal t,
                                                  const PetscReal x[],
                                                  PetscInt numComponents,
                                                  PetscScalar *s,
                                                  void *context)
    {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(1 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = trace_strain(x[0], x[1], t);

        return 0;
    } // solnkernel_trace_strain

    static PetscErrorCode solnkernel_displacement_t(PetscInt spaceDim,
                                                    PetscReal t,
                                                    const PetscReal x[],
                                                    PetscInt numComponents,
                                                    PetscScalar *s,
                                                    void *context)
    {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = displacement_t_x(x[0], x[1], t);
        s[1] = displacement_t_y(x[0], x[1], t);

        return 0;
    } // solnkernel_displacement_t

    static PetscErrorCode solnkernel_pressure_t(PetscInt spaceDim,
                                                PetscReal t,
                                                const PetscReal x[],
                                                PetscInt numComponents,
                                                PetscScalar *s,
                                                void *context)
    {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(1 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = pressure_t(x[0], x[1], t);

        return 0;
    } // solnkernel_pressure_t

    static PetscErrorCode solnkernel_trace_strain_t(PetscInt spaceDim,
                                                    PetscReal t,
                                                    const PetscReal x[],
                                                    PetscInt numComponents,
                                                    PetscScalar *s,
                                                    void *context)
    {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(1 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = trace_strain_t(x[0], x[1], t);

        return 0;
    } // solnkernel_trace_strain_t

protected:
    void setUp(void)
    {
        TestIsotropicLinearPoroelasticity::setUp();

        // Overwrite component names for control of debugging info at test level.
        GenericComponent::setName("TestIsotropicLinearPoroelasticity2D_TS_LT");
        pythia::journal::debug_t debug(GenericComponent::getName());
        // debug.activate(); // DEBUGGING

        CPPUNIT_ASSERT(!_data);
        _data = new TestPoroelasticity_Data();
        CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->spaceDim = 2;
        _data->meshFilename = ":UNKNOWN:"; // Set in child class.
        _data->boundaryLabel = "boundary";

        CPPUNIT_ASSERT(!_data->cs);
        _data->cs = new spatialdata::geocoords::CSCart;
        CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->spaceDim);

        CPPUNIT_ASSERT(_data->normalizer);
        _data->normalizer->setLengthScale(LENGTHSCALE);
        _data->normalizer->setTimeScale(TIMESCALE);
        _data->normalizer->setPressureScale(PRESSURESCALE);
        _data->normalizer->computeDensityScale();

        _data->startTime = 0.0;
        _data->endTime = 0.1;
        _data->timeStep = 0.05;

        // solnDiscretizations set in derived class.

        _data->numAuxSubfields = 9;
        static const char *_auxSubfields[9] = {
            // order must match order of subfields in auxiliary field
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
        _data->auxSubfields = _auxSubfields;
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
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_auxDiscretizations);

        CPPUNIT_ASSERT(_data->auxDB);
        _data->auxDB->addValue("solid_density", solid_density, solid_density_units());
        _data->auxDB->addValue("fluid_density", fluid_density, fluid_density_units());
        _data->auxDB->addValue("fluid_viscosity", fluid_viscosity, fluid_viscosity_units());
        _data->auxDB->addValue("porosity", porosity, porosity_units());
        _data->auxDB->addValue("shear_modulus", shear_modulus, shear_modulus_units());
        _data->auxDB->addValue("drained_bulk_modulus", drained_bulk_modulus, drained_bulk_modulus_units());
        _data->auxDB->addValue("biot_coefficient", biot_coefficient, biot_coefficient_units());
        _data->auxDB->addValue("fluid_bulk_modulus", fluid_bulk_modulus, fluid_bulk_modulus_units());
        _data->auxDB->addValue("solid_bulk_modulus", solid_bulk_modulus, solid_bulk_modulus_units());
        _data->auxDB->addValue("isotropic_permeability", isotropic_permeability, isotropic_permeability_units());
        _data->auxDB->setCoordSys(*_data->cs);

        CPPUNIT_ASSERT(_material);
        _material->setFormulation(pylith::problems::Physics::QUASISTATIC);
        _material->useBodyForce(false);
        _material->useSourceDensity(false);
        _rheology->useReferenceState(false);

        _rheology->useTensorPermeability(false);

        _material->setDescriptiveLabel("Isotropic Linear Poroelasticity Plane Strain");
        _material->setMaterialId(24);

        // Displacement BC
        static const PylithInt constrainedDisplacementDOF[2] = {0, 1};
        static const PylithInt numConstrainedDisplacement = 2;
        _bcDisplacement->setConstrainedDOF(constrainedDisplacementDOF, numConstrainedDisplacement);
        _bcDisplacement->setMarkerLabel("boundary");
        _bcDisplacement->setSubfieldName("displacement");
        _bcDisplacement->setUserFn(solnkernel_displacement);
        _bcDisplacement->setUserFnDot(solnkernel_displacement_t);

        // Pressure BC
        static const PylithInt constrainedPressureDOF[1] = {0};
        static const PylithInt numConstrainedPressure = 1;
        _bcPressure->setConstrainedDOF(constrainedPressureDOF, numConstrainedPressure);
        _bcPressure->setMarkerLabel("boundary");
        _bcPressure->setSubfieldName("pressure");
        _bcPressure->setUserFn(solnkernel_pressure);
        _bcPressure->setUserFnDot(solnkernel_pressure_t);

    } // setUp

    // Set exact solution in domain.
    void _setExactSolution(void)
    {
        CPPUNIT_ASSERT(_solution);
        PetscErrorCode err = 0;
        PetscDS prob = NULL;
        PetscWeakForm wf = NULL;
        DMLabel label = NULL;
        err = DMGetDS(_solution->getDM(), &prob);
        CPPUNIT_ASSERT(!err);
        err = DMGetLabel(_solution->getDM(), "material-id", &label);
        CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_displacement, NULL);
        CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolutionTimeDerivative(prob, 0, solnkernel_displacement_t, NULL);
        CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 1, solnkernel_pressure, NULL);
        CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolutionTimeDerivative(prob, 1, solnkernel_pressure_t, NULL);
        CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 2, solnkernel_trace_strain, NULL);
        CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolutionTimeDerivative(prob, 2, solnkernel_trace_strain_t, NULL);
        CPPUNIT_ASSERT(!err);
        err = PetscDSGetWeakForm(prob, &wf);
        CPPUNIT_ASSERT(!err);
        err = PetscWeakFormSetIndexResidual(wf, label, 24, 0, 0, 1, pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0_mms_tl_u, 1, NULL);
        CPPUNIT_ASSERT(!err);
        err = PetscWeakFormSetIndexResidual(wf, label, 24, 1, 0, 1, pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0_mms_tl_p, 1, NULL);
        CPPUNIT_ASSERT(!err);
    } // _setExactSolution

}; // TestIsotropicLinearPoroelasticity2D_TS_LT
const double pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT::LENGTHSCALE = 1.0;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT::TIMESCALE = 1.0;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT::PRESSURESCALE = 1.0;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT::BODYFORCE = 1.0;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT::XMIN = -4.0;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT::XMAX = 4.0;

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT_TriP2 : public pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT
{
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity2D_TS_LT_TriP2, TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void)
    {
        TestIsotropicLinearPoroelasticity2D_TS_LT::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 3;
        static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
            pylith::topology::Field::Discretization(2, 2), // displacement
            pylith::topology::Field::Discretization(1, 2), // pressure
            pylith::topology::Field::Discretization(1, 2), // trace_strain
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_solnDiscretizations);

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
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity2D_TS_LT_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT_TriP3 : public pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT
{
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity2D_TS_LT_TriP3,
                           TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void)
    {
        TestIsotropicLinearPoroelasticity2D_TS_LT::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 3;
        static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
            pylith::topology::Field::Discretization(3, 3), // displacement
            pylith::topology::Field::Discretization(2, 3), // pressure
            pylith::topology::Field::Discretization(2, 3), // trace_strain
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_solnDiscretizations);

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
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity2D_TS_LT_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT_TriP4 : public pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT
{
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity2D_TS_LT_TriP4, TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void)
    {
        TestIsotropicLinearPoroelasticity2D_TS_LT::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 3;
        static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
            pylith::topology::Field::Discretization(4, 4), // displacement
            pylith::topology::Field::Discretization(3, 4), // pressure
            pylith::topology::Field::Discretization(3, 4), // trace_strain
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
            pylith::topology::Field::Discretization(0, 4), // solid_density
            pylith::topology::Field::Discretization(0, 4), // fluid_density
            pylith::topology::Field::Discretization(0, 4), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 4), // porosity
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // drained_bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // biot_coefficient
            pylith::topology::Field::Discretization(0, 4), // biot_modulus
            pylith::topology::Field::Discretization(0, 4), // isotropic_permeability
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity2D_TS_LT_TriP4
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT_TriP4);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ2 : public pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT
{
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ2, TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void)
    {
        TestIsotropicLinearPoroelasticity2D_TS_LT::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        _data->numSolnSubfields = 3;
        static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
            pylith::topology::Field::Discretization(2, 2), // displacement
            pylith::topology::Field::Discretization(1, 2), // pressure
            pylith::topology::Field::Discretization(1, 2), // trace_strain
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_solnDiscretizations);

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
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ3 : public pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT
{
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ3, TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void)
    {
        TestIsotropicLinearPoroelasticity2D_TS_LT::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        _data->numSolnSubfields = 3;
        static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
            pylith::topology::Field::Discretization(3, 3), // displacement
            pylith::topology::Field::Discretization(2, 3), // pressure
            pylith::topology::Field::Discretization(2, 3), // trace_strain
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_solnDiscretizations);

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
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ4 : public pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT
{
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ4, TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void)
    {
        TestIsotropicLinearPoroelasticity2D_TS_LT::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        _data->numSolnSubfields = 3;
        static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
            pylith::topology::Field::Discretization(4, 4), // displacement
            pylith::topology::Field::Discretization(3, 4), // pressure
            pylith::topology::Field::Discretization(3, 4), // trace_strain
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
            pylith::topology::Field::Discretization(0, 4), // solid_density
            pylith::topology::Field::Discretization(0, 4), // fluid_density
            pylith::topology::Field::Discretization(0, 4), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 4), // porosity
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // drained_bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // biot_coefficient
            pylith::topology::Field::Discretization(0, 4), // biot_modulus
            pylith::topology::Field::Discretization(0, 4), // isotropic_permeability
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ4
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity2D_TS_LT_QuadQ4);

// End of file
