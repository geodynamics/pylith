// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "tests/src/MMSTest.hh" // ISA MMSTEST

#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity
#include "pylith/sources/MomentTensorForce.hh" // USES MomentTensorForce
#include "pylith/sources/RickerWavelet.hh" // USES RickerWavelet
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "pylith/scales/Scales.hh" // USES Scales

#include "pylith/problems/Physics.hh" // USES FormulationEnum
#include "pylith/topology/Field.hh" // HASA FieldBase::Discretization

namespace pylith {
    class TestMomentTensorSource;
    class TestMomentTensorSource_Data;
}

/// C++ class for testing moment tensor source using MMS.
class pylith::TestMomentTensorSource : public pylith::testing::MMSTest {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] data Data for MMS test.
     */
    TestMomentTensorSource(TestMomentTensorSource_Data* data);

    /// Destructor.
    ~TestMomentTensorSource(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Initialize objects for test.
    void _initialize(void);

    /// Set exact solution and time derivative of solution in domain.
    void _setExactSolution(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    TestMomentTensorSource_Data* _data; ///< Test parameters.

}; // class TestMomentTensorSource

// ================================================================================================
class pylith::TestMomentTensorSource_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestMomentTensorSource_Data(void);

    /// Destructor
    ~TestMomentTensorSource_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    const char* journalName; ///< Name for MMSTest journals.
    int spaceDim; ///< Spatial dimension of domain.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    const char* meshOptions; ///< Command line options for mesh.
    const char* boundaryLabel; ///< Group defining domain boundary.
    bool useAsciiMesh; ///< Use MeshIOAscii to read mesh, otherwise use PETSc.

    PylithReal jacobianConvergenceRate; ///< Expected convergence rate for Jacobian (when not linear).
    PylithReal tolerance; ///< Tolerance for discretization and residual test.
    bool isJacobianLinear; ///< Jacobian should be linear.
    bool allowZeroResidual; ///< Allow residual to be exactly zero.

    PylithReal t; ///< Time for MMS solution.
    PylithReal dt; ///< Time step in simulation.
    spatialdata::geocoords::CSCart cs; ///< Coordinate system.
    pylith::scales::Scales scales; ///< Scales for nondimensionalization.
    pylith::problems::Physics::FormulationEnum formulation; ///< Time stepping formulation

    pylith::materials::Elasticity material; ///< Material.
    pylith::materials::IsotropicLinearElasticity rheology; ///< Bulk rheology for material.
    pylith::sources::MomentTensorForce source; ///< Moment tensor source.
    pylith::sources::RickerWavelet sourceTimeFunction; ///< Source time function.
    std::vector<pylith::bc::BoundaryCondition*> bcs; ///< Dirichlet boundary conditions.

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

    // Source auxiliary fields.
    size_t numSourceAuxSubfields; ///< Number of auxiliary subfields for source.
    const char** sourceAuxSubfields; ///< Names of auxiliary subfields for source.
    pylith::topology::Field::Discretization const* sourceAuxDiscretizations; ///< Discretizations for source auxiliary subfields.
    spatialdata::spatialdb::UserFunctionDB sourceAuxDB; ///< Spatial database for source auxiliary field.

}; // TestMomentTensorSource_Data

// End of file
