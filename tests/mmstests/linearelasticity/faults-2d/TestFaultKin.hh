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

/**
 * @file mmstets/linearelasticity/faults-2d/TestFaultKin.hh
 *
 * @brief C++ class for testing faults with prescribed slip.
 */

#if !defined(pylith_mmstests_testfaultkin_hh)
#define pylith_mmstests_testfaultkin_hh

#include "tests/src/MMSTest.hh" // ISA MMSTEST

#include "pylith/faults/faultsfwd.hh" // HOLDSA FaultCohesiveKin
#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/problems/Physics.hh" // USES FormulationEnum
#include "pylith/topology/Field.hh" // HASA FieldBase::Discretization

namespace pylith {
    class TestFaultKin;
    class TestFaultKin_Data;
} // pylith

class pylith::TestFaultKin : public pylith::testing::MMSTest {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] data Data for MMS test.
     */
    TestFaultKin(TestFaultKin_Data* data);

    /// Destructor.
    ~TestFaultKin(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Initialize objects for test.
    void _initialize(void);

    /// Set exact solution and time derivative of solution in domain.
    void _setExactSolution(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    TestFaultKin_Data* _data; ///< Test parameters.

}; // class TestFaultKin

// ================================================================================================
class pylith::TestFaultKin_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestFaultKin_Data(void);

    /// Destructor
    ~TestFaultKin_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    const char* journalName; ///< Name for MMSTest journals.
    int spaceDim; ///< Spatial dimension of domain.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    const char* meshOptions; ///< Command line options for mesh.
    const char* boundaryLabel; ///< Group defining domain boundary.
    bool useAsciiMesh; ///< Use MeshIOAscii to read mesh, otherwise use PETSc.

    PylithReal jacobianConvergenceRate; ///< Expected convergence rate for Jacobiab (when not linear).
    PylithReal tolerance; ///< Tolerance for discretization and residual test.
    bool isJacobianLinear; ///< Jacobian is should be linear.
    bool allowZeroResidual; ///< Allow residual to be exactly zero.

    PylithReal t; ///< Time for MMS solution.
    PylithReal dt; ///< Time step in simulation.
    spatialdata::geocoords::CSCart cs; ///< Coordinate system.
    spatialdata::units::Nondimensional normalizer; ///< Scales for nondimensionalization.
    pylith::problems::Physics::FormulationEnum formulation; ///< Time stepping formulation

    std::vector<pylith::materials::Material*> materials; ///< Materials.
    pylith::materials::IsotropicLinearElasticity rheology; ///< Bulk rheology for materials.
    std::vector<pylith::bc::BoundaryCondition*> bcs; ///< Dirichlet boundary condition.
    std::vector<pylith::faults::FaultCohesive*> faults; ///< Fault interface conditions.
    spatialdata::spatialdb::GravityField* gravityField; ///< Gravity field.

    // Solution field.
    int numSolnSubfieldsDomain; ///< Number of solution fields for domain.
    int numSolnSubfieldsFault; ///< Number of solution fields for fault.
    pylith::topology::Field::Discretization const* solnDiscretizations; ///< Discretizations for solution fields.

    /// Array of functions providing exact solution.
    pylith::testing::MMSTest::solution_fn* exactSolnFns;

    /// Array of functions providing exact solution time derivative.
    pylith::testing::MMSTest::solution_fn* exactSolnDotFns;

    // Material auxiliary fields.
    int matNumAuxSubfields;
    const char** matAuxSubfields;
    pylith::topology::Field::Discretization const* matAuxDiscretizations;
    spatialdata::spatialdb::UserFunctionDB matAuxDB;

    // Fault auxiliary fields.
    int faultNumAuxSubfields;
    const char** faultAuxSubfields;
    pylith::topology::Field::Discretization const* faultAuxDiscretizations;
    spatialdata::spatialdb::UserFunctionDB faultAuxDB;

    pylith::faults::KinSrc* kinSrc;

}; // TestFaultKin_Data

#endif // pylith_mmstests_testfaultkin_hh

// End of file
