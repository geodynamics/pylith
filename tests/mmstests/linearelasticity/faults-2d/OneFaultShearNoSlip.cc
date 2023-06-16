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

/** @file tests/mmstests/faults/TestFaultKin2D_OneFaultShearNoSlip.cc
 *
 * Square domain of sides 12.0 km with a single through-going fault at y=+2 km.
 * Simple shear with no slip.
 */

#include <portinfo>

#include "OneFaultShearNoSlip.hh" // Implementation of cases

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/KinSrcStep.hh" // USES KinSrcStep
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn
#include "pylith/bc/NeumannUserFn.hh" // USES NeumannUserFn

#include "pylith/topology/Mesh.hh" // USES pylith::topology::Mesh::cells_label_name
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

namespace pylith {
    class _OneFaultShearNoSlip;
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::_OneFaultShearNoSlip {
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

        return strain_xy() * 2.0 * mu / 2.25e+10;
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
        r0[1] += tractionShear*tanDir[1] + tractionNormal*n[1];
    } // boundary_tractions

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

    static PetscErrorCode solnkernel_lagrangemultiplier(PetscInt spaceDim,
                                                        PetscReal t,
                                                        const PetscReal x[],
                                                        PetscInt numComponents,
                                                        PetscScalar* s,
                                                        void* context) {
        assert(2 == spaceDim);
        assert(x);
        assert(2 == numComponents);
        assert(s);

        s[0] = faulttraction_x(x[0], x[1]);
        s[1] = faulttraction_y(x[0], x[1]);

        return PETSC_SUCCESS;
    } // solnkernel_lagrangemultiplier

public:

    static
    TestFaultKin_Data* createData(void) {
        TestFaultKin_Data* data = new TestFaultKin_Data();assert(data);

        data->journalName = "OneFaultShearNoSLip";

        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.

        data->normalizer.setLengthScale(1.0e+03);
        data->normalizer.setTimeScale(2.0);
        data->normalizer.setPressureScale(2.25e+10);
        data->normalizer.computeDensityScale();

        // solnDiscretizations set in derived class.

        data->matNumAuxSubfields = 3;
        static const char* _matAuxSubfields[3] = {"density", "shear_modulus", "bulk_modulus"};
        data->matAuxSubfields = _matAuxSubfields;
        static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        data->matAuxDB.addValue("density", density, density_units());
        data->matAuxDB.addValue("vp", vp, vp_units());
        data->matAuxDB.addValue("vs", vs, vs_units());
        data->matAuxDB.setCoordSys(data->cs);

        assert(!data->kinSrc);
        data->kinSrc = new pylith::faults::KinSrcStep();assert(data->kinSrc);
        data->kinSrc->setOriginTime(0.0);
        data->faultAuxDB.addValue("initiation_time", initiation_time, time_units());
        data->faultAuxDB.addValue("final_slip_opening", finalslip_opening, slip_units());
        data->faultAuxDB.addValue("final_slip_left_lateral", finalslip_leftlateral, slip_units());
        data->faultAuxDB.setCoordSys(data->cs);

        data->faultNumAuxSubfields = 1;
        static const char* _faultAuxSubfields[1] = { "slip" };
        data->faultAuxSubfields = _faultAuxSubfields;
        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 1), // slip
        };
        data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        // Materials
        data->materials.resize(3);
        { // xneg
            pylith::materials::Elasticity* material = new pylith::materials::Elasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setDescription("Isotropic Linear Elasticity Plane Strain");
            material->setLabelValue(10);
            material->setBulkRheology(&data->rheology);
            data->materials[0] = material;
        } // xneg
        { // mid
            pylith::materials::Elasticity* material = new pylith::materials::Elasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setDescription("Isotropic Linear Elasticity Plane Strain");
            material->setLabelValue(20);
            material->setBulkRheology(&data->rheology);
            data->materials[1] = material;
        } // mid
        { // xpos
            pylith::materials::Elasticity* material = new pylith::materials::Elasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setDescription("Isotropic Linear Elasticity Plane Strain");
            material->setLabelValue(15);
            material->setBulkRheology(&data->rheology);
            data->materials[2] = material;
        } // xpos

        // Boundary conditions
        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        data->bcs.resize(4);
        { // boundary_xneg
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[0] = bc;
        } // boundary_xneg
        { // boundary_xpos
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[1] = bc;
        } // boundary_xpos
        { // boundary_yneg
            pylith::bc::NeumannUserFn* bc = new pylith::bc::NeumannUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_yneg");
            bc->setLabelValue(1);
            bc->setUserFn(boundary_tractions);
            data->bcs[2] = bc;
        } // boundary_yneg
        { // boundary_ypos
            pylith::bc::NeumannUserFn* bc = new pylith::bc::NeumannUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_ypos");
            bc->setLabelValue(1);
            bc->setUserFn(boundary_tractions);
            data->bcs[3] = bc;
        } // boundary_ypos

        // Faults
        data->faults.resize(1);
        { // xpos
            pylith::faults::FaultCohesiveKin* fault = new pylith::faults::FaultCohesiveKin();
            fault->setCohesiveLabelValue(100);
            fault->setSurfaceLabelName("fault_xpos");

            const int numRuptures = 1;
            const char* ruptureNames[1] = { "rupture" };
            pylith::faults::KinSrc* ruptures[1] = { data->kinSrc };
            fault->setEqRuptures(ruptureNames, numRuptures, ruptures, numRuptures);
            data->faults[0] = fault;
        } // xpos

        pylith::utils::PetscOptions options;
        options.add("-fieldsplit_displacement_pc_type", "lu");
        options.override ();

        data->numSolnSubfieldsDomain = 1;
        data->numSolnSubfieldsFault = 1;
        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[2] = {
            solnkernel_disp,
            solnkernel_lagrangemultiplier,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        data->exactSolnDotFns = nullptr;

        return data;
    } // createData

}; // TestFaultKin2D_OneFaultShearNoSlip

// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::OneFaultShearNoSlip::TriP1(void) {
    TestFaultKin_Data* data = pylith::_OneFaultShearNoSlip::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    assert(1 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP1


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::OneFaultShearNoSlip::TriP2(void) {
    TestFaultKin_Data* data = pylith::_OneFaultShearNoSlip::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 2), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    assert(1 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP2


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::OneFaultShearNoSlip::TriP3(void) {
    TestFaultKin_Data* data = pylith::_OneFaultShearNoSlip::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 3), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    assert(1 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(3, 3, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP3


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::OneFaultShearNoSlip::TriP4(void) {
    TestFaultKin_Data* data = pylith::_OneFaultShearNoSlip::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 4), // density
        pylith::topology::Field::Discretization(0, 4), // shear_modulus
        pylith::topology::Field::Discretization(0, 4), // bulk_modulus
    };
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 4), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    assert(1 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(4, 4), // disp
        pylith::topology::Field::Discretization(4, 4, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP4


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::OneFaultShearNoSlip::QuadQ1(void) {
    TestFaultKin_Data* data = pylith::_OneFaultShearNoSlip::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    assert(1 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ1


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::OneFaultShearNoSlip::QuadQ2(void) {
    TestFaultKin_Data* data = pylith::_OneFaultShearNoSlip::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 2), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    assert(1 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ2


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::OneFaultShearNoSlip::QuadQ3(void) {
    TestFaultKin_Data* data = pylith::_OneFaultShearNoSlip::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 3), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    assert(1 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(3, 3, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ3


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::OneFaultShearNoSlip::QuadQ4(void) {
    TestFaultKin_Data* data = pylith::_OneFaultShearNoSlip::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 4), // density
        pylith::topology::Field::Discretization(0, 4), // shear_modulus
        pylith::topology::Field::Discretization(0, 4), // bulk_modulus
    };
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 4), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    assert(1 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(4, 4), // disp
        pylith::topology::Field::Discretization(4, 4, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ4


// End of file
