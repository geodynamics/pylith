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
// Copyright (c) 2010-2023 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file tests/mmstests/linearelasticity/faults-2d/UniformStrainRateExplicit.cc
 *
 * The domain is two cells, one on each side of a fault. The solution corresponds to a constant
 * left-lateral slip rate of 1.5 m/s with fixed boundaries.
 */

#include <portinfo>

#include "UniformStrainRateExplicit.hh" // ISA UniformStrainRateExplicit

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

// #define USE_QUASISTATIC

// ------------------------------------------------------------------------------------------------
namespace pylith {
    class _UniformStrainRateExplicit;
} // pylith

class pylith::_UniformStrainRateExplicit {
    // Solution parameters
    // Strain rates
    static const double ExxDotN;
    static const double EyyDotN;
    static const double ExyDotN;
    static const double ExxDotP;
    static const double EyyDotP;
    static const double ExyDotP;
    static const double SLIP_RATE_N; // nondimensional
    static const double TIME_SNAPSHOT_N; // nondimensional
    static const double X_FAULT_N; // nondimensional

    // Dimensionless
    static const double LENGTH_SCALE;
    static const double TIME_SCALE;
    static const double PRESSURE_SCALE;

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
        const double vsN = 3000.0;
        return (x < X_FAULT_N) ? vsN : (fabs(ExyDotP) > 0.0) ? vsN * sqrt(ExyDotN / ExyDotP) : vsN;
    } // vs

    static const char* vs_units(void) {
        return "m/s";
    } // vs_units

    // Vp
    static double vp(const double x,
                     const double y) {
        const double vpN = sqrt(3.0) * vs(X_FAULT_N-1.0, y);
        double value = vpN;
        if (x > X_FAULT_N) {
            const double densityN = density(X_FAULT_N-1.0, y);
            const double densityP = density(X_FAULT_N+1.0, y);

            const double vsN = vs(X_FAULT_N-1.0, y);
            const double muN = vsN*vsN * densityN;

            const double vsP = vs(X_FAULT_N+1.0, y);
            const double muP = vsP*vsP * densityP;

            const double lambdaN = vpN*vpN * densityN - 2.0*muN;
            const double lambdaP = (fabs(ExxDotP + EyyDotP) > 0.0) ? (lambdaN*(ExxDotN + EyyDotN) + 2.0*muN*ExxDotN - 2.0*muP*ExxDotP) /  (ExxDotP + EyyDotP) : lambdaN;
            const double vpP = sqrt((lambdaP + 2.0*muP) / densityP);
            value = vpP;
        } // if
        return value;
    } // vp

    static const char* vp_units(void) {
        return "m/s";
    } // vp_units

    static double lameMuNondim(const double x,
                               const double y) {
        const double vsXY = vs(x, y);
        const double densityXY = density(x, y);
        const double mu = vsXY*vsXY * densityXY / PRESSURE_SCALE;
        return mu;
    }

    static double lameLambdaNondim(const double x,
                                   const double y) {
        const double mu = lameMuNondim(x, y);
        const double vpXY = vp(x, y);
        const double densityXY = density(x, y);
        const double lambda = (vpXY*vpXY * densityXY) / PRESSURE_SCALE - 2.0*mu;
        return lambda;
    }

    // Kinematic rupture auxiliary components.

    // Initiation time
    static double initiation_time(const double x,
                                  const double y) {
        return 0.0 * TIME_SCALE;
    } // initiation_time

    static const char* time_units(void) {
        return "s";
    } // time_units

    // Slip rate
    static double finalslip_opening(const double x,
                                    const double y) {
        return 0.0;
    } // finalslip_opening

    static double finalslip_leftlateral(const double x,
                                        const double y) {
        return SLIP_RATE_N * TIME_SNAPSHOT_N * LENGTH_SCALE;
    } // finalslip_leftlateral

    static const char* slip_units(void) {
        return "m";
    } // slip_units

    // Solution subfields (nondimensional).

    // Acceleration
    static double acc_x(const double x,
                        const double y) {
        return 0.0;
    } // acc_x

    static double acc_y(const double x,
                        const double y) {
        return 0.0;
    } // acc_y

    static const char* acc_units(void) {
        return "m/s*2";
    } // vel_units

    // Velocity
    static double vel_x(const double x,
                        const double y,
                        PetscInt flag) {
        double vel = 0;
        if ((flag < 0) || (x < X_FAULT_N)) {
            vel = ExxDotN * (x-X_FAULT_N) + 0*ExyDotN * y;
        } else {
            vel = ExxDotP * (x-X_FAULT_N) + 0*ExyDotP * y;
        }
        return vel;
    } // vel_x

    static double vel_y(const double x,
                        const double y,
                        PetscInt flag) {
        double vel = 0.0;

        if ((flag < 0) || (x < X_FAULT_N)) {
            vel = 2.0*ExyDotN * (x-X_FAULT_N) + EyyDotN * y + 0.5 * SLIP_RATE_N;
        } else {
            vel = 2.0*ExyDotP * (x-X_FAULT_N) + EyyDotP * y - 0.5 * SLIP_RATE_N;
        } // if/else
        return vel;
    } // vel_y

    static const char* vel_units(void) {
        return "m/s";
    } // vel_units

    // Displacement
    static double disp_x(const double x,
                         const double y,
                         PetscInt flag) {
        return vel_x(x, y, flag) * TIME_SNAPSHOT_N;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         PetscInt flag) {
        const double disp = vel_y(x, y, flag) * TIME_SNAPSHOT_N;
        return disp;
    } // disp_y

    static double faulttraction_x(const double x,
                                  const double y) {
        const double muN = lameMuNondim(X_FAULT_N-1.0, y);
        const double lambdaN = lameLambdaNondim(X_FAULT_N-1.0, y);
        const double stress_xx = (lambdaN * (ExxDotN + EyyDotN) + 2.0*muN * ExxDotN) * TIME_SNAPSHOT_N;
        return stress_xx;
    } // faulttraction_x

    static double faulttraction_y(const double x,
                                  const double y) {
        const double muN = lameMuNondim(X_FAULT_N-1.0, y);
        const double stress_xy = 2.0 * muN * ExyDotN * TIME_SNAPSHOT_N;
        return stress_xy;
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

        const double muN = lameMuNondim(X_FAULT_N-1.0, x[1]);

        const PylithScalar tanDir[2] = {-n[1], n[0] };
        const double strain_xy = ExyDotN * TIME_SNAPSHOT_N;
        const PylithScalar tractionShear = -strain_xy * 2.0 * muN;
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
                PYLITH_JOURNAL_LOGICERROR("Unknown cell type in solution displacement kernel.");
            }
            flag = cell < numCellsLeftFault ? -1 : +1;
        } // if
        s[0] = disp_x(x[0], x[1], flag);
        s[1] = disp_y(x[0], x[1], flag);

        return PETSC_SUCCESS;
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
                PYLITH_JOURNAL_LOGICERROR("Unknown cell type in solution displacement kernel.");
            }
            flag = cell < numCellsLeftFault ? -1 : +1;
        } // if
        s[0] = vel_x(x[0], x[1], flag);
        s[1] = vel_y(x[0], x[1], flag);

        return PETSC_SUCCESS;
    } // solnkernel_disp

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

        s[0] = acc_x(x[0], x[1]);
        s[1] = acc_y(x[0], x[1]);

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

        data->journalName = "UniformStrainRateExplicit";
        data->allowZeroResidual = true;
        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.

        data->normalizer.setLengthScale(LENGTH_SCALE);
        data->normalizer.setTimeScale(TIME_SCALE);
        data->normalizer.setPressureScale(2.25e+10);
        data->normalizer.computeDensityScale();

#if !defined(USE_QUASISTATIC)
        data->formulation = pylith::problems::Physics::DYNAMIC_IMEX;
#else
        data->formulation = pylith::problems::Physics::QUASISTATIC;
#endif
        data->t = TIME_SNAPSHOT_N;
        data->dt = 0.05;

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
            material->setFormulation(data->formulation);
            material->useBodyForce(false);
            material->setDescription("Isotropic Linear Elasticity Plane Strain");
            material->setLabelValue(10);
            material->setBulkRheology(&data->rheology);
            data->materials[0] = material;
        } // xneg
        { // mid
            pylith::materials::Elasticity* material = new pylith::materials::Elasticity();assert(material);
            material->setFormulation(data->formulation);
            material->useBodyForce(false);
            material->setDescription("Isotropic Linear Elasticity Plane Strain");
            material->setLabelValue(20);
            material->setBulkRheology(&data->rheology);
            data->materials[1] = material;
        } // mid
        { // xpos
            pylith::materials::Elasticity* material = new pylith::materials::Elasticity();assert(material);
            material->setFormulation(data->formulation);
            material->useBodyForce(false);
            material->setDescription("Isotropic Linear Elasticity Plane Strain");
            material->setLabelValue(15);
            material->setBulkRheology(&data->rheology);
            data->materials[2] = material;
        } // xpos

        // Boundary conditions
        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
#if !defined(USE_QUASISTATIC)
        data->bcs.resize(6);
        { // boundary_xneg displacement
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(solnkernel_disp);
            bc->setUserFnDot(solnkernel_vel);
            data->bcs[0] = bc;
        } // boundary_xneg displacement
        { // boundary_xneg velocity
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("velocity");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(solnkernel_vel);
            bc->setUserFnDot(solnkernel_acc);
            data->bcs[1] = bc;
        } // boundary_xneg velocity
        { // boundary_xpos displacement
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(solnkernel_disp);
            bc->setUserFnDot(solnkernel_vel);
            data->bcs[2] = bc;
        } // boundary_xpos displacement
        { // boundary_xpos velocity
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("velocity");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(solnkernel_vel);
            bc->setUserFnDot(solnkernel_acc);
            data->bcs[3] = bc;
        } // boundary_xpos velocity
        { // boundary_yneg traction
            pylith::bc::NeumannUserFn* bc = new pylith::bc::NeumannUserFn();
            bc->setSubfieldName("velocity");
            bc->setLabelName("boundary_yneg");
            bc->setLabelValue(1);
            bc->setUserFn(boundary_tractions);
            data->bcs[4] = bc;
        } // boundary_yneg traction
        { // boundary_ypos traction
            pylith::bc::NeumannUserFn* bc = new pylith::bc::NeumannUserFn();
            bc->setSubfieldName("velocity");
            bc->setLabelName("boundary_ypos");
            bc->setLabelValue(1);
            bc->setUserFn(boundary_tractions);
            data->bcs[5] = bc;
        } // boundary_ypos traction
#else
        data->bcs.resize(4);
        { // boundary_xneg displacement
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[0] = bc;
        } // boundary_xneg displacement
        { // boundary_xpos displacement
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[1] = bc;
        } // boundary_xpos displacement
        { // boundary_yneg traction
            pylith::bc::NeumannUserFn* bc = new pylith::bc::NeumannUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_yneg");
            bc->setLabelValue(1);
            bc->setUserFn(boundary_tractions);
            data->bcs[2] = bc;
        } // boundary_yneg traction
        { // boundary_ypos traction
            pylith::bc::NeumannUserFn* bc = new pylith::bc::NeumannUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_ypos");
            bc->setLabelValue(1);
            bc->setUserFn(boundary_tractions);
            data->bcs[3] = bc;
        } // boundary_ypos traction
#endif
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

#if !defined(USE_QUASISTATIC)
        data->numSolnSubfieldsDomain = 2;
        data->numSolnSubfieldsFault = 1;
        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[3] = {
            solnkernel_disp,
            solnkernel_vel,
            solnkernel_lagrangemultiplier,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        static const pylith::testing::MMSTest::solution_fn _exactSolnDotFns[3] = {
            solnkernel_vel,
            solnkernel_acc,
            solnkernel_lagrangemultiplier,
        };
        data->exactSolnDotFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnDotFns);
#else
        data->numSolnSubfieldsDomain = 1;
        data->numSolnSubfieldsFault = 1;
        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[2] = {
            solnkernel_disp,
            solnkernel_lagrangemultiplier,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        data->exactSolnDotFns = nullptr;
#endif

        return data;
    } // createData

}; // TestFaultKin2D_UniformStrainRateExplicit

const double pylith::_UniformStrainRateExplicit::ExxDotN = 0.0;
const double pylith::_UniformStrainRateExplicit::EyyDotN = 0.0;
const double pylith::_UniformStrainRateExplicit::ExyDotN = -0.1;
const double pylith::_UniformStrainRateExplicit::ExxDotP = 0.0;
const double pylith::_UniformStrainRateExplicit::EyyDotP = 0.0;
const double pylith::_UniformStrainRateExplicit::ExyDotP = -0.1;
const double pylith::_UniformStrainRateExplicit::SLIP_RATE_N = 0.0;// 3.0e-3;
const double pylith::_UniformStrainRateExplicit::TIME_SNAPSHOT_N = 5.0;
const double pylith::_UniformStrainRateExplicit::X_FAULT_N = +2.0;

const double pylith::_UniformStrainRateExplicit::LENGTH_SCALE = 1000.0;
const double pylith::_UniformStrainRateExplicit::PRESSURE_SCALE = 2.25e+10;
const double pylith::_UniformStrainRateExplicit::TIME_SCALE = 2.0;

// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::UniformStrainRateExplicit::TriP1(void) {
    TestFaultKin_Data* data = pylith::_UniformStrainRateExplicit::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // vel
        pylith::topology::Field::Discretization(1, 1, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP1


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::UniformStrainRateExplicit::TriP2(void) {
    TestFaultKin_Data* data = pylith::_UniformStrainRateExplicit::createData();assert(data);

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

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2), // vel
        pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP2


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::UniformStrainRateExplicit::TriP3(void) {
    TestFaultKin_Data* data = pylith::_UniformStrainRateExplicit::createData();assert(data);

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

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(3, 3), // vel
        pylith::topology::Field::Discretization(3, 3, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP3


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::UniformStrainRateExplicit::TriP4(void) {
    TestFaultKin_Data* data = pylith::_UniformStrainRateExplicit::createData();assert(data);

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

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(4, 4), // disp
        pylith::topology::Field::Discretization(4, 4), // vel
        pylith::topology::Field::Discretization(4, 4, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP4


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::UniformStrainRateExplicit::QuadQ1(void) {
    TestFaultKin_Data* data = pylith::_UniformStrainRateExplicit::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

#if !defined(USE_QUASISTATIC)
    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // vel
        pylith::topology::Field::Discretization(1, 1, 1, -1, true), // lagrange_multiplier_fault
    };
#else
    assert(1 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1, 1, -1, true), // lagrange_multiplier_fault
    };
#endif
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ1


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::UniformStrainRateExplicit::QuadQ2(void) {
    TestFaultKin_Data* data = pylith::_UniformStrainRateExplicit::createData();assert(data);

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

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2), // vel
        pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ2


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::UniformStrainRateExplicit::QuadQ3(void) {
    TestFaultKin_Data* data = pylith::_UniformStrainRateExplicit::createData();assert(data);

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

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(3, 3), // vel
        pylith::topology::Field::Discretization(3, 3, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ3


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKin_Data*
pylith::UniformStrainRateExplicit::QuadQ4(void) {
    TestFaultKin_Data* data = pylith::_UniformStrainRateExplicit::createData();assert(data);

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

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(4, 4), // disp
        pylith::topology::Field::Discretization(4, 4), // vel
        pylith::topology::Field::Discretization(4, 4, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ4


// End of file
