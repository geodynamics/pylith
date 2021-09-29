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

/** @file tests/mmstests/faults/TestFaultKin2D_ShearNoSlipStatic.cc
 *
 * Square domain of sides 8.0 km with a through-going fault running
 * through the center in the y-direction. Simple shear with no slip.
 */

#include <portinfo>

#include "TestFaultKin.hh" // ISA TestFaultKin

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/KinSrcStep.hh" // USES KinSrcStep
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn
#include "pylith/bc/NeumannUserFn.hh" // USES NeumannUserFn

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

namespace pylith {
    namespace mmstests {
        class TestFaultKin2D_ShearNoSlipStatic;

        class TestFaultKin2D_ShearNoSlipStatic_TriP1;
        class TestFaultKin2D_ShearNoSlipStatic_TriP2;
        class TestFaultKin2D_ShearNoSlipStatic_TriP3;
        class TestFaultKin2D_ShearNoSlipStatic_TriP4;

        class TestFaultKin2D_ShearNoSlipStatic_QuadQ1;
        class TestFaultKin2D_ShearNoSlipStatic_QuadQ2;
        class TestFaultKin2D_ShearNoSlipStatic_QuadQ3;
        class TestFaultKin2D_ShearNoSlipStatic_QuadQ4;
    } // tests/mmstests
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic :
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
        return 0.0;
    } // slip_leftlateral

    static const char* slip_units(void) {
        return "m";
    } // slip_units

    // Solution subfields.
    static double strain_xx(void) {
        return 0.0;
    } // strain_xx

    static double strain_yy(void) {
        return 0.0;
    } // strain_yy

    static double strain_xy(void) {
        return 0.3;
    } // strain_xy

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return strain_xx()*x + strain_xy()*y;
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return strain_xy()*x + strain_yy()*y;
    } // disp_y

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    static double faulttraction_x(const double x,
                                  const double y) {
        return 0.0;
    } // faulttraction_x

    static double faulttraction_y(const double x,
                                  const double y) {
        const double mu = density(x, y) * vs(x, y) * vs(x, y);

        return -strain_xy() * 2.0 * mu / 2.25e+10;
    } // faulttraction_y

    static
    void boundary_tractions(const PylithInt dim,
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
                            const PylithReal x[],
                            const PylithReal n[],
                            const PylithInt numConstants,
                            const PylithScalar constants[],
                            PylithScalar r0[]) {
        assert(r0);
        const double mu = density(x[0], x[1]) * vs(x[0], x[1]) * vs(x[0], x[1]);

        const PylithScalar tanDir[2] = {-n[1], n[0] };
        const PylithScalar tractionShear = -strain_xy() * 2.0 * mu / 2.25e+10;
        const PylithScalar tractionNormal = 0.0;
        r0[0] += tractionShear*tanDir[0] + tractionNormal*n[0];
        r0[0] += tractionShear*tanDir[1] + tractionNormal*n[1];
    } // boundary_tractions

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
        GenericComponent::setName("TestFaultKin2D_ShearNoSlipStatic");

        CPPUNIT_ASSERT(!_data);
        _data = new TestFaultKin_Data();CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->spaceDim = 2;
        _data->meshFilename = ":UNKNOWN:"; // Set in child class.

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

        // Materials
        _materials.resize(2);
        { // xneg
            pylith::materials::Elasticity* material = new pylith::materials::Elasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setDescriptiveLabel("Isotropic Linear Elasticity Plane Strain");
            material->setMaterialId(10);
            material->setBulkRheology(_data->rheology);
            _materials[0] = material;
        } // xneg
        { // xpos
            pylith::materials::Elasticity* material = new pylith::materials::Elasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setDescriptiveLabel("Isotropic Linear Elasticity Plane Strain");
            material->setMaterialId(20);
            material->setBulkRheology(_data->rheology);
            _materials[1] = material;
        } // xpos

        // Boundary conditions
        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        _bcs.resize(4);
        { // boundary_xneg
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setMarkerLabel("boundary_xneg");
            bc->setSubfieldName("displacement");
            bc->setUserFn(solnkernel_disp);
            _bcs[0] = bc;
        } // boundary_xneg
        { // boundary_xpos
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setMarkerLabel("boundary_xpos");
            bc->setSubfieldName("displacement");
            bc->setUserFn(solnkernel_disp);
            _bcs[1] = bc;
        } // boundary_xpos
        { // boundary_yneg
            pylith::bc::NeumannUserFn* bc = new pylith::bc::NeumannUserFn();
            bc->setMarkerLabel("boundary_yneg");
            bc->setSubfieldName("displacement");
            bc->setUserFn(boundary_tractions);
            _bcs[2] = bc;
        } // boundary_yneg
        { // boundary_ypos
            pylith::bc::NeumannUserFn* bc = new pylith::bc::NeumannUserFn();
            bc->setMarkerLabel("boundary_ypos");
            bc->setSubfieldName("displacement");
            bc->setUserFn(boundary_tractions);
            _bcs[3] = bc;
        } // boundary_ypos

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

}; // TestFaultKin2D_ShearNoSlipStatic

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_TriP1 :
    public pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ShearNoSlipStatic_TriP1,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ShearNoSlipStatic::setUp();
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

}; // TestFaultKin2D_ShearNoSlipStatic_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_TriP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_TriP2 :
    public pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ShearNoSlipStatic_TriP2,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ShearNoSlipStatic::setUp();
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

}; // TestFaultKin2D_ShearNoSlipStatic_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_TriP3 :
    public pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ShearNoSlipStatic_TriP3,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ShearNoSlipStatic::setUp();
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

}; // TestFaultKin2D_ShearNoSlipStatic_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_TriP4 :
    public pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ShearNoSlipStatic_TriP4,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ShearNoSlipStatic::setUp();
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

}; // TestFaultKin2D_ShearNoSlipStatic_TriP4
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_TriP4);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_QuadQ1 :
    public pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ShearNoSlipStatic_QuadQ1,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ShearNoSlipStatic::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(1, 1), // disp
            pylith::topology::Field::Discretization(1, 1, 1, -1, true), // lagrange_multiplier_fault
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKin2D_ShearNoSlipStatic_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_QuadQ1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_QuadQ2 :
    public pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ShearNoSlipStatic_QuadQ2,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ShearNoSlipStatic::setUp();
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

}; // TestFaultKin2D_ShearNoSlipStatic_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_QuadQ3 :
    public pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ShearNoSlipStatic_QuadQ3,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ShearNoSlipStatic::setUp();
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

}; // TestFaultKin2D_ShearNoSlipStatic_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_QuadQ3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_QuadQ4 :
    public pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKin2D_ShearNoSlipStatic_QuadQ4,
                           TestFaultKin);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKin2D_ShearNoSlipStatic::setUp();
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

}; // TestFaultKin2D_ShearNoSlipStatic_QuadQ4
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKin2D_ShearNoSlipStatic_QuadQ4);

// End of file
