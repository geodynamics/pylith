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
 * @file mmstets/linearporoelasticity/nofaults-2d/TestLinearPoroElasticity.hh
 *
 * @brief C++ class for testing PoroElasticity for various rheologies.
 */

#if !defined(pylith_mmstests_testlinearporoelasticity_hh)
#define pylith_mmstests_testlinearporoelasticity_hh

#include "tests/src/MMSTest.hh" // ISA MMSTEST

#include "pylith/materials/Poroelasticity.hh" // USES Poroelasticity
#include "pylith/materials/IsotropicLinearPoroelasticity.hh" // USES IsotropicLinearPoroelasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/problems/Physics.hh" // USES FormulationEnum
#include "pylith/topology/Field.hh" // HASA FieldBase::Discretization

namespace pylith {
    class TestLinearPoroElasticity;
    class TestLinearPoroElasticity_Data;
}

/// C++ class for testing PoroElasticity using a variety of MMS tests.
class pylith::TestLinearPoroElasticity : public pylith::testing::MMSTest {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] data Data for MMS test.
     */
    TestLinearPoroElasticity(TestLinearPoroElasticity_Data* data);

    /// Destructor.
    ~TestLinearPoroElasticity(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Initialize objects for test.
    void _initialize(void);

    /// Set exact solution and time derivative of solution in domain.
    void _setExactSolution(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    TestLinearPoroElasticity_Data* _data; ///< Test parameters.

}; // class TestLinearElasticity

// ================================================================================================
class pylith::TestLinearPoroElasticity_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestLinearPoroElasticity_Data(void);

    /// Destructor
    ~TestLinearPoroElasticity_Data(void);

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

    pylith::materials::Poroelasticity material; ///< Materials.
    pylith::materials::IsotropicLinearPoroelasticity rheology; ///< Bulk rheology for materials.
    spatialdata::spatialdb::GravityField* gravityField; ///< Gravity field.
    std::vector<pylith::bc::BoundaryCondition*> bcs; ///< Dirichlet boundary condition.

    // Solution field.
    size_t numSolnSubfields; ///< Number of solution fields.
    pylith::topology::Field::Discretization const* solnDiscretizations; ///< Discretizations for solution fields.

    /// Array of functions providing exact solution.
    pylith::testing::MMSTest::solution_fn* exactSolnFns;

    /// Array of functions providing exact solution time derivative.
    pylith::testing::MMSTest::solution_fn* exactSolnDotFns;

    // Material auxiliary fields.
    size_t numAuxSubfields; ///< Number of auxiliary subfields for materials.
    const char** auxSubfields; ///< Names of auxiliary subfields for materials.
    pylith::topology::Field::Discretization const* auxDiscretizations; ///< Discretizations for auxiliary subfields.
    spatialdata::spatialdb::UserFunctionDB auxDB; ///< Spatial database for auxiliary field.

}; // TestLinearPoroElasticity_Data

#endif // pylith_testlinearporoelasticity_hh

// End of file
