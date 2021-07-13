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

/** @file tests/mmstests/faults/TestFaultKin2D_RigidBlocksStatic.cc
 *
 * Square domain of sides 8.0 km with a through-going fault running
 * through the center in the y-direction. The two opposing sides each
 * move as rigid blocks with 3.0 m of right-lateral slip.
 */

#include <portinfo>

#include "TestFaultKin.hh" // ISA TestFaultKin

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/KinSrcStep.hh" // USES KinSrcStep
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

namespace pylith {
    namespace mmstests {
        class TestFaultKin2D_RigidBlocksStatic;

        class TestFaultKin2D_RigidBlocksStatic_TriP1;
        class TestFaultKin2D_RigidBlocksStatic_TriP2;
        class TestFaultKin2D_RigidBlocksStatic_TriP3;
        class TestFaultKin2D_RigidBlocksStatic_TriP4;

        class TestFaultKin2D_RigidBlocksStatic_QuadQ1;
        class TestFaultKin2D_RigidBlocksStatic_QuadQ2;
        class TestFaultKin2D_RigidBlocksStatic_QuadQ3;
        class TestFaultKin2D_RigidBlocksStatic_QuadQ4;
    } // tests/mmstests
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_RigidBlocksStatic :
    public pylith::mmstests::TestFaultKin {
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
        return sqrt(3.0)*vs(x,y);
    } // vp

    static const char* vp_units(void) {
        return "m/s";
    } // vp_units

    // Kinematic rupture auxiliary components.

    // Initiation time
    static double initiation_time(const double x,
                                  const double y) {
        return 0.0;
    } // initiation_time

    static const char* time_units(void) {
        return "s";
    } // time_units

    // Slip
    static double finalslip_opening(const double x,
                                    const double y) {
        return 0.0;
    } // slip_opening

    static double finalslip_leftlateral(const double x,
                                        const double y) {
        return -1.5;
    } // slip_leftlateral

    static const char* slip_units(void) {
        return "m";
    } // slip_units

    // Solution subfields.

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return 0.0;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         PetscInt flag) {
        const double disp = 0.75e-3;
        if (!flag) {
            return x < 0.0 ? +disp : -disp;
        } else {
            return flag < 0 ? +disp : -disp;
        } // if/else
    } // disp_y

    static double faulttraction_x(const double x,
                                  const double y) {
        return 0.0;
    } // faulttraction_x

    static double faulttraction_y(const double x,
                                  const double y) {
        return 0.0;
    } // faulttraction_y

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
        PetscInt flag = 0;
        if (context) {
            PetscInt cell = 0;
            DMPolytopeType cellType = DM_POLYTOPE_UNKNOWN;
            DMPlexGetActivePoint((PetscDM) context, &cell);
            DMPlexGetCellType((PetscDM) context, cell, &cellType);
            PetscInt numCellsLeftFault = 0;
            switch (cellType) {
            case DM_POLYTOPE_TRIANGLE:
                numCellsLeftFault = 4;
                break;
            case DM_POLYTOPE_QUADRILATERAL:
                numCellsLeftFault = 2;
                break;
            default:
                CPPUNIT_FAIL("Unknown cell type in solution displacement kernel.");
            }
            flag = cell < numCellsLeftFault ? -1 : +1;
        } // if
        s[1] = disp_y(x[0], x[1], flag);

        return 0;
    } // solnkernel_disp

    static PetscErrorCode solnkernel_lagrangemultiplier(PetscInt spaceDim,
                                                        PetscReal t,
                                                        const PetscReal x[],
                                                        PetscInt numComponents,
                                                        PetscScalar* s,
                                                        void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = faulttraction_x(x[0], x[1]);
        s[1] = faulttraction_y(x[0], x[1]);

        return 0;
    } // solnkernel_lagrangemultiplier

protected:

    void setUp(void) {
        TestFaultKin::setUp();

        // Overwrite component name for control of journals at test level.
        GenericComponent::setName("TestFaultKin2D_RigidBlocksStatic");

        CPPUNIT_ASSERT(!_data);
        _data = new TestFaultKin_Data();CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

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

        _data->matNumAuxSubfields = 3;
        static const char* _matAuxSubfields[3] = {"density", "shear_modulus", "bulk_modulus"};
        _data->matAuxSubfields = _matAuxSubfields;
        static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        CPPUNIT_ASSERT(!_data->rheology);
        _data->rheology = new pylith::materials::IsotropicLinearElasticity();CPPUNIT_ASSERT(_data->rheology);
        CPPUNIT_ASSERT(_data->matAuxDB);
        _data->matAuxDB->addValue("density", density, density_units());
        _data->matAuxDB->addValue("vp", vp, vp_units());
        _data->matAuxDB->addValue("vs", vs, vs_units());
        _data->matAuxDB->setCoordSys(*_data->cs);

        CPPUNIT_ASSERT(!_data->kinsrc);
        _data->kinsrc = new pylith::faults::KinSrcStep();CPPUNIT_ASSERT(_data->kinsrc);
        _data->kinsrc->originTime(0.0);
        CPPUNIT_ASSERT(_data->faultAuxDB);
        _data->faultAuxDB->addValue("initiation_time", initiation_time, time_units());
        _data->faultAuxDB->addValue("final_slip_opening", finalslip_opening, slip_units());
        _data->faultAuxDB->addValue("final_slip_left_lateral", finalslip_leftlateral, slip_units());
        _data->faultAuxDB->setCoordSys(*_data->cs);

        _data->faultNumAuxSubfields = 1;
        static const char* _faultAuxSubfields[1] = { "slip" };
        _data->faultAuxSubfields = _faultAuxSubfields;
        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 1), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        CPPUNIT_ASSERT(_material);
        _material->setFormulation(pylith::problems::Physics::QUASISTATIC);
        _material->useBodyForce(false);
        _material->setDescriptiveLabel("Isotropic Linear Elascitity Plane Strain");
        _material->setMaterialId(24);
        _material->setBulkRheology(_data->rheology);

        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        _bc->setConstrainedDOF(constrainedDOF, numConstrained);
        _bc->setMarkerLabel("boundary");
        _bc->setSubfieldName("displacement");
        _bc->setUserFn(solnkernel_disp);

        _fault->setInterfaceId(100);
        _fault->setSurfaceMarkerLabel("fault");

    } // setUp

    // Set exact solution in domain.
    void _setExactSolution(void) {
        CPPUNIT_ASSERT(_solution);

        PetscDM dm = _solution->getDM();
        PetscDMLabel label;
        PetscIS is;
        PetscInt cohesiveCell;
        PetscErrorCode err = 0;
        PetscDS prob = NULL;
        err = DMGetDS(dm, &prob);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_disp, dm);CPPUNIT_ASSERT(!err);
        err = DMGetLabel(dm, "material-id", &label);CPPUNIT_ASSERT(!err);
        err = DMLabelGetStratumIS(label, _fault->getInterfaceId(), &is);CPPUNIT_ASSERT(!err);
        err = ISGetMinMax(is, &cohesiveCell, NULL);CPPUNIT_ASSERT(!err);
        err = ISDestroy(&is);CPPUNIT_ASSERT(!err);
        err = DMGetCellDS(dm, cohesiveCell, &prob);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_disp, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 1, solnkernel_lagrangemultiplier, NULL);CPPUNIT_ASSERT(!err);
    } // _setExactSolution

}; // TestFaultKin2D_RigidBlocksStatic

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_TriP1 :
    public pylith::mmstests::TestFaultKin2D_RigidBlocksStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_RigidBlocksStatic_TriP1,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_RigidBlocksStatic::setUp();
        CPPUNIT_ASSERT(_data);
        _allowZeroResidual = true;

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(1, 1), // disp
            pylith::topology::Field::Discretization(1, 1, 1, -1, true), // lagrange_multiplier_fault
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKin2D_RigidBlocksStatic_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_TriP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_TriP2 :
    public pylith::mmstests::TestFaultKin2D_RigidBlocksStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_RigidBlocksStatic_TriP2,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_RigidBlocksStatic::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 2), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKin2D_RigidBlocksStatic_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_TriP3 :
    public pylith::mmstests::TestFaultKin2D_RigidBlocksStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_RigidBlocksStatic_TriP3,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_RigidBlocksStatic::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 3), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(3, 3), // disp
            pylith::topology::Field::Discretization(3, 3, 1, -1, true), // lagrange_multiplier_fault
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKin2D_RigidBlocksStatic_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_TriP4 :
    public pylith::mmstests::TestFaultKin2D_RigidBlocksStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_RigidBlocksStatic_TriP4,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_RigidBlocksStatic::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 4), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(4, 4), // disp
            pylith::topology::Field::Discretization(4, 4, 1, -1, true), // lagrange_multiplier_fault
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKin2D_RigidBlocksStatic_TriP4
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_TriP4);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_QuadQ1 :
    public pylith::mmstests::TestFaultKin2D_RigidBlocksStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_RigidBlocksStatic_QuadQ1,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_RigidBlocksStatic::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(1, 1), // disp
            pylith::topology::Field::Discretization(1, 1, 1, -1, true), // lagrange_multiplier_fault
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKin2D_RigidBlocksStatic_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_QuadQ1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_QuadQ2 :
    public pylith::mmstests::TestFaultKin2D_RigidBlocksStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_RigidBlocksStatic_QuadQ2,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_RigidBlocksStatic::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 2), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKin2D_RigidBlocksStatic_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_QuadQ3 :
    public pylith::mmstests::TestFaultKin2D_RigidBlocksStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_RigidBlocksStatic_QuadQ3,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_RigidBlocksStatic::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 3), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(3, 3), // disp
            pylith::topology::Field::Discretization(3, 3, 1, -1, true), // lagrange_multiplier_fault
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKin2D_RigidBlocksStatic_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_QuadQ3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_QuadQ4 :
    public pylith::mmstests::TestFaultKin2D_RigidBlocksStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_RigidBlocksStatic_QuadQ4,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_RigidBlocksStatic::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 4), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(4, 4), // disp
            pylith::topology::Field::Discretization(4, 4, 1, -1, true), // lagrange_multiplier_fault
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKin2D_RigidBlocksStatic_QuadQ4
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_RigidBlocksStatic_QuadQ4);

// End of file
