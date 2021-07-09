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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestIsotropicLinearIncompElasticity.hh" // ISA TestIsotropicLinearIncompElasticity2D

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/IncompressibleElasticity.hh" // USES IncompressibleElasticity
#include "pylith/materials/IsotropicLinearIncompElasticity.hh" // USES IsotropicLinearIncompElasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

namespace pylith {
    namespace mmstests {
        class TestIsotropicLinearIncompElasticity2D_UniformPressure;

        class TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP1;
        class TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP2;
        class TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP3;
        class TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP4;

        class TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ1;
        class TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ2;
        class TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ3;
        class TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ4;
    } // tests/mmstests
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure :
    public pylith::mmstests::TestIsotropicLinearIncompElasticity {
    // Spatial database user functions for auxiiliary subfields (includes derived fields).

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

    // Solution subfields.

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return 0.0;
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return 0.0;
    } // disp_y

    // Pressure
    static double pressure(const double x,
                           const double y) {
        return 5.0+6;
    } // pressure

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

        s[0] = disp_x(x[0], x[1]);
        s[1] = disp_y(x[0], x[1]);

        return 0;
    } // solnkernel_disp

    static PetscErrorCode solnkernel_pressure(PetscInt spaceDim,
                                              PetscReal t,
                                              const PetscReal x[],
                                              PetscInt numComponents,
                                              PetscScalar* s,
                                              void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(1 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = pressure(x[0], x[1]);

        return 0;
    } // solnkernel_pressure

protected:

    void setUp(void) {
        TestIsotropicLinearIncompElasticity::setUp();

        // Overwrite component names for control of journals at test level.
        GenericComponent::setName("TestIsotropicLinearIncompElasticity2D_UniformPressure");

        CPPUNIT_ASSERT(!_data);
        _data = new TestIncompressibleElasticity_Data();CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;
        _jacobianConvergenceRate = 1.0;

        _data->spaceDim = 2;
        _data->meshFilename = ":UNKNOWN:"; // Set in child class.
        _data->boundaryLabel = "boundary";

        CPPUNIT_ASSERT(!_data->cs);
        _data->cs = new spatialdata::geocoords::CSCart;CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->spaceDim);

        CPPUNIT_ASSERT(_data->normalizer);
        _data->normalizer->setLengthScale(1.0e+03);
        _data->normalizer->setTimeScale(2.0);
        _data->normalizer->setPressureScale(2.25e+10);
        _data->normalizer->computeDensityScale();

        _data->startTime = 0.0;
        _data->endTime = 0.1;
        _data->timeStep = 0.05;

        // solnDiscretizations set in derived class.

        _data->numAuxSubfields = 3;
        static const char* _auxSubfields[3] = {"density", "shear_modulus", "bulk_modulus"};
        _data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        CPPUNIT_ASSERT(_data->auxDB);
        _data->auxDB->addValue("density", density, density_units());
        _data->auxDB->addValue("vp", vp, vp_units());
        _data->auxDB->addValue("vs", vs, vs_units());
        _data->auxDB->setCoordSys(*_data->cs);

        CPPUNIT_ASSERT(_material);
        _material->useBodyForce(false);
        _rheology->useReferenceState(false);

        _material->setDescriptiveLabel("Isotropic Linear Incompressible Elascitity Plane Strain");
        _material->setMaterialId(24);

        static const PylithInt constrainedDispDOF[2] = {0, 1};
        static const PylithInt numConstrainedDisp = 2;
        _bcDisplacement->setConstrainedDOF(constrainedDispDOF, numConstrainedDisp);
        _bcDisplacement->setMarkerLabel("boundary");
        _bcDisplacement->setSubfieldName("displacement");
        _bcDisplacement->setUserFn(solnkernel_disp);

        static const PylithInt constrainedPressureDOF[1] = {0};
        static const PylithInt numConstrainedPressure = 1;
        _bcPressure->setConstrainedDOF(constrainedPressureDOF, numConstrainedPressure);
        _bcPressure->setMarkerLabel("boundary");
        _bcPressure->setSubfieldName("pressure");
        _bcPressure->setUserFn(solnkernel_pressure);

    } // setUp

    // Set exact solution in domain.
    void _setExactSolution(void) {
        CPPUNIT_ASSERT(_solution);

        PetscErrorCode err = 0;
        PetscDS prob = NULL;
        err = DMGetDS(_solution->getDM(), &prob);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_disp, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 1, solnkernel_pressure, NULL);CPPUNIT_ASSERT(!err);
    } // _setExactSolution

}; // TestIsotropicLinearIncompElasticity2D_UniformPressure

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP1 :
    public pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP1,
                           TestIsotropicLinearIncompElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearIncompElasticity2D_UniformPressure::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(1, 1), // disp
            pylith::topology::Field::Discretization(0, 1), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP2 :
    public pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP2,
                           TestIsotropicLinearIncompElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearIncompElasticity2D_UniformPressure::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(1, 2), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP3 :
    public pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP3,
                           TestIsotropicLinearIncompElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearIncompElasticity2D_UniformPressure::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(3, 3), // disp
            pylith::topology::Field::Discretization(2, 3), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP4 :
    public pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP4,
                           TestIsotropicLinearIncompElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearIncompElasticity2D_UniformPressure::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(4, 4), // disp
            pylith::topology::Field::Discretization(3, 4), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP4
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_TriP4);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ1 :
    public pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ1,
                           TestIsotropicLinearIncompElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearIncompElasticity2D_UniformPressure::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(1, 1), // disp
            pylith::topology::Field::Discretization(0, 1), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ2 :
    public pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ2,
                           TestIsotropicLinearIncompElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearIncompElasticity2D_UniformPressure::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(1, 2), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ3 :
    public pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ3,
                           TestIsotropicLinearIncompElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearIncompElasticity2D_UniformPressure::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(3, 3), // disp
            pylith::topology::Field::Discretization(2, 3), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ4 :
    public pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ4,
                           TestIsotropicLinearIncompElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearIncompElasticity2D_UniformPressure::setUp();
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
            pylith::topology::Field::Discretization(3, 4), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ4
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearIncompElasticity2D_UniformPressure_QuadQ4);

// End of file
