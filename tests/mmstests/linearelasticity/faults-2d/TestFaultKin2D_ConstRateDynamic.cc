// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file tests/mmstests/faults/TestFaultKin2D_ConstRateDynamic.cc
 *
 * The domain is two cells, one on each side of a fault. The solution corresponds to a constant
 * left-lateral slip rate of 1.5 m/s with fixed boundaries.
 */

#include <portinfo>

#include "TestFaultKin.hh" // ISA TestFaultKin

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/KinSrcConstRate.hh" // USES KinSrcConstRate
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "pylith/topology/Mesh.hh" // USES pylith::topology::Mesh::cells_label_name
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Scales.hh" // USES Scales

namespace pylith {
    namespace mmstests {
        class TestFaultKin2D_ConstRateDynamic;

        class TestFaultKin2D_ConstRateDynamic_TriP1;
        class TestFaultKin2D_ConstRateDynamic_TriP2;
        class TestFaultKin2D_ConstRateDynamic_TriP3;
        class TestFaultKin2D_ConstRateDynamic_TriP4;

        class TestFaultKin2D_ConstRateDynamic_QuadQ1;
        class TestFaultKin2D_ConstRateDynamic_QuadQ2;
        class TestFaultKin2D_ConstRateDynamic_QuadQ3;
        class TestFaultKin2D_ConstRateDynamic_QuadQ4;
    } // tests/mmstests
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ConstRateDynamic :
    public pylith::mmstests::TestFaultKin {
    // Spatial database user functions for auxiiliary subfields (includes derived fields).

    // Dimensionless
    static const double SLIPRATE;
    static const double VELOCITY;
    static const double TIMESTAMP;
    static const double LENGTH_SCALE;
    static const double TIME_SCALE;
    static const double PRESSURE_SCALE;
    static const double DOMAIN_X;

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

    // Slip rate
    static double sliprate_opening(const double x,
                                   const double y) {
        return 0.0;
    } // sliprate_opening

    static double sliprate_leftlateral(const double x,
                                       const double y) {
        std::cout << "slip rate="<<SLIPRATE * LENGTH_SCALE/TIME_SCALE<<std::endl;

        return SLIPRATE * LENGTH_SCALE/TIME_SCALE;
    } // sliprate_leftlateral

    static const char* sliprate_units(void) {
        return "m/s";
    } // sliprate_units

    // Solution subfields (nondimensional).

    // Velocity
    static double velocity_x(const double x,
                             const double y) {
        return 0.0;
    } // velocity_x

    static double velocity_y(const double x,
                             const double y,
                             PetscInt flag) {
        double amplitude = 0.0;
        if (!flag) {
            amplitude = x < 0 ? -0.5*SLIPRATE : +0.5*SLIPRATE;
        } else {
            amplitude = flag < 0 ? -0.5*SLIPRATE : +0.5*SLIPRATE;
        } // if/else
        amplitude += (VELOCITY - 0.5*SLIPRATE) * x / (0.5*DOMAIN_X);
        return amplitude;
    } // velocity_y

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return 0.0;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         PetscInt flag) {
        const double disp = velocity_y(x, y, flag) * TIMESTAMP;
        return disp;
    } // disp_y

    static double faulttraction_x(const double x,
                                  const double y) {
        return 0.0;
    } // faulttraction_x

    static double faulttraction_y(const double x,
                                  const double y) {
        const double shearModulusN = density(x,y) * pow(vs(x,y), 2) / PRESSURE_SCALE;

        return shearModulusN * (VELOCITY - 0.5*SLIPRATE) * TIMESTAMP / (0.5*DOMAIN_X);
    } // faulttraction_y

    static PetscErrorCode bckernel_disp(PetscInt spaceDim,
                                        PetscReal t,
                                        const PetscReal x[],
                                        PetscInt numComponents,
                                        PetscScalar* s,
                                        void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = 0.0;
        s[1] = disp_y(x[0], x[1], 0);
        std::cout << "BOUNDARY x=("<<x[0]<<", "<<x[1]<<"), s=("<<s[0]<<", "<<s[1]<<")" << std::endl;

        return 0;
    } // bckernel_disp

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
                numCellsLeftFault = 12;
                break;
            case DM_POLYTOPE_QUADRILATERAL:
                numCellsLeftFault = 6;
                break;
            default:
                CPPUNIT_FAIL("Unknown cell type in solution displacement kernel.");
            }
            flag = cell < numCellsLeftFault ? -1 : +1;
        } // if
        s[1] = disp_y(x[0], x[1], flag);
        std::cout << "SOLUTION x=("<<x[0]<<", "<<x[1]<<"), s=("<<s[0]<<", "<<s[1]<<")" << std::endl;

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
        std::cout << "x=("<<x[0]<<", "<<x[1]<<"), traction=("<<faulttraction_x(x[0], x[1])<<", "<<faulttraction_y(x[0], x[1])<<")" << std::endl;

        s[0] = faulttraction_x(x[0], x[1]);
        s[1] = faulttraction_y(x[0], x[1]);

        return 0;
    } // solnkernel_lagrangemultiplier

protected:

    void setUp(void) {
        TestFaultKin::setUp();

        // Overwrite component name for control of journals at test level.
        GenericComponent::setName("TestFaultKin2D_ConstRateDynamic");

        CPPUNIT_ASSERT(!_data);
        _data = new TestFaultKin_Data();CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->spaceDim = 2;
        _data->meshFilename = ":UNKNOWN:"; // Set in child class.

        CPPUNIT_ASSERT(!_data->cs);
        _data->cs = new spatialdata::geocoords::CSCart;CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->spaceDim);

        CPPUNIT_ASSERT(_data->scales);
        _data->scales->setLengthScale(LENGTH_SCALE);
        _data->scales->setTimeScale(TIME_SCALE);
        _data->scales->setPressureScale(2.25e+10);

        _data->startTime = 0.0;
        _data->endTime = 2.0*TIMESTAMP;
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
        _data->kinsrc = new pylith::faults::KinSrcConstRate();CPPUNIT_ASSERT(_data->kinsrc);
        CPPUNIT_ASSERT(_data->faultAuxDB);
        _data->faultAuxDB->addValue("initiation_time", initiation_time, time_units());
        _data->faultAuxDB->addValue("slip_rate_opening", sliprate_opening, sliprate_units());
        _data->faultAuxDB->addValue("slip_rate_left_lateral", sliprate_leftlateral, sliprate_units());
        _data->faultAuxDB->setCoordSys(*_data->cs);

        _data->faultNumAuxSubfields = 2;
        static const char* _faultAuxSubfields[2] = { "slip", "slip acceleration" };
        _data->faultAuxSubfields = _faultAuxSubfields;
        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[2] = {
            pylith::topology::Field::Discretization(0, 1), // slip
            pylith::topology::Field::Discretization(0, 1), // slip acceleration
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        // Materials
        _materials.resize(3);
        { // xneg
            pylith::materials::Elasticity* material = new pylith::materials::Elasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setIdentifier("elasticity");
            material->setName("material-id=10");
            material->setLabelValue(10);
            material->setBulkRheology(_data->rheology);
            _materials[0] = material;
        } // xneg
        { // mid
            pylith::materials::Elasticity* material = new pylith::materials::Elasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setIdentifier("elasticity");
            material->setName("material-id=20");
            material->setLabelValue(20);
            material->setBulkRheology(_data->rheology);
            _materials[1] = material;
        } // mid
        { // xpos
            pylith::materials::Elasticity* material = new pylith::materials::Elasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setIdentifier("elasticity");
            material->setName("material-id=15");
            material->setLabelValue(15);
            material->setBulkRheology(_data->rheology);
            _materials[2] = material;
        } // xpos

        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        _bcs.resize(2);
        { // boundary_xpos
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(bckernel_disp);
            _bcs[0] = bc;
        } // boundary_xpos
        { // boundary_xneg
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(bckernel_disp);
            _bcs[1] = bc;
        } // boundary_zneg

        // Faults
        _faults.resize(1);
        { // xpos
            pylith::faults::FaultCohesiveKin* fault = new pylith::faults::FaultCohesiveKin();
            fault->setCohesiveLabelValue(100);
            fault->setSurfaceLabelName("fault_xpos");

            const int numRuptures = 1;
            const char* ruptureNames[1] = { "rupture" };
            pylith::faults::KinSrc* ruptures[1] = { _data->kinsrc };
            fault->setEqRuptures(ruptureNames, numRuptures, ruptures, numRuptures);
            _faults[0] = fault;
        } // xpos

    } // setUp

    // Set exact solution in domain.
    void _setExactSolution(void) {
        const pylith::topology::Field* solution = _problem->getSolution();
        CPPUNIT_ASSERT(solution);

        PetscDM dm = solution->getDM();
        PetscDMLabel label;
        PetscIS is;
        PetscInt cohesiveCell;
        PetscErrorCode err = 0;
        PetscDS ds = NULL;
        err = DMGetDS(dm, &ds);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(ds, 0, solnkernel_disp, dm);CPPUNIT_ASSERT(!err);
        err = DMGetLabel(dm, pylith::topology::Mesh::cells_label_name, &label);CPPUNIT_ASSERT(!err);
        err = DMLabelGetStratumIS(label, _faults[0]->getCohesiveLabelValue(), &is);CPPUNIT_ASSERT(!err);
        err = ISGetMinMax(is, &cohesiveCell, NULL);CPPUNIT_ASSERT(!err);
        err = ISDestroy(&is);CPPUNIT_ASSERT(!err);
        err = DMGetCellDS(dm, cohesiveCell, &ds, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(ds, 0, solnkernel_disp, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(ds, 1, solnkernel_lagrangemultiplier, NULL);CPPUNIT_ASSERT(!err);
    } // _setExactSolution

}; // TestFaultKin2D_ConstRateDynamic
#if 0
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::SLIPRATE = 1.5;
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::VELOCITY = 4.0;
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::TIMESTAMP = 10.0;
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::LENGTH_SCALE = 1000.0;
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::PRESSURE_SCALE = 2.5e+10;
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::TIME_SCALE = 2.0;
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::DOMAIN_X = 8000.0;
#else
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::SLIPRATE = 3.0;
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::VELOCITY = 1.5;
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::TIMESTAMP = 10.0;
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::LENGTH_SCALE = 1.0;
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::PRESSURE_SCALE = 2.0e+6;
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::TIME_SCALE = 2.0;
const double pylith::mmstests::TestFaultKin2D_ConstRateDynamic::DOMAIN_X = 8000.0;

#endif

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ConstRateDynamic_TriP1 :
    public pylith::mmstests::TestFaultKin2D_ConstRateDynamic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ConstRateDynamic_TriP1,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ConstRateDynamic::setUp();
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

}; // TestFaultKin2D_ConstRateDynamic_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ConstRateDynamic_TriP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ConstRateDynamic_TriP2 :
    public pylith::mmstests::TestFaultKin2D_ConstRateDynamic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ConstRateDynamic_TriP2,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ConstRateDynamic::setUp();
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

}; // TestFaultKin2D_ConstRateDynamic_TriP2
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ConstRateDynamic_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ConstRateDynamic_TriP3 :
    public pylith::mmstests::TestFaultKin2D_ConstRateDynamic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ConstRateDynamic_TriP3,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ConstRateDynamic::setUp();
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

}; // TestFaultKin2D_ConstRateDynamic_TriP3
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ConstRateDynamic_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ConstRateDynamic_TriP4 :
    public pylith::mmstests::TestFaultKin2D_ConstRateDynamic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ConstRateDynamic_TriP4,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ConstRateDynamic::setUp();
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

}; // TestFaultKin2D_ConstRateDynamic_TriP4
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ConstRateDynamic_TriP4);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ConstRateDynamic_QuadQ1 :
    public pylith::mmstests::TestFaultKin2D_ConstRateDynamic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ConstRateDynamic_QuadQ1,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ConstRateDynamic::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(1, 1), // disp
            pylith::topology::Field::Discretization(1, 1, 1, -1, true), // lagrange_multiplier_fault
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKin2D_ConstRateDynamic_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ConstRateDynamic_QuadQ1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ConstRateDynamic_QuadQ2 :
    public pylith::mmstests::TestFaultKin2D_ConstRateDynamic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ConstRateDynamic_QuadQ2,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ConstRateDynamic::setUp();
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

}; // TestFaultKin2D_ConstRateDynamic_QuadQ2
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ConstRateDynamic_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ConstRateDynamic_QuadQ3 :
    public pylith::mmstests::TestFaultKin2D_ConstRateDynamic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ConstRateDynamic_QuadQ3,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ConstRateDynamic::setUp();
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

}; // TestFaultKin2D_ConstRateDynamic_QuadQ3
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ConstRateDynamic_QuadQ3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ConstRateDynamic_QuadQ4 :
    public pylith::mmstests::TestFaultKin2D_ConstRateDynamic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ConstRateDynamic_QuadQ4,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ConstRateDynamic::setUp();
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

}; // TestFaultKin2D_ConstRateDynamic_QuadQ4
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ConstRateDynamic_QuadQ4);

// End of file
