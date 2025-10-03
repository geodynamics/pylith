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

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/problems/problemsfwd.hh" // HOLDSA SolutionFactory
#include "pylith/topology/Field.hh" // HOLDSA Field::SubfieldInfo
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA Coordsys
#include "pylith/scales/scalesfwd.hh" // HOLDSA Scales

#include <map> // USES std::map

/// Namespace for pylith package
namespace pylith {
    namespace problems {
        class TestSolutionFactory;
        class TestSolutionFactory_Data;
    } // problems
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::problems::TestSolutionFactory : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestSolutionFactory(TestSolutionFactory_Data* data);

    /// Destructor.
    ~TestSolutionFactory(void);

    /// Test adding displacement and velocity subfields.
    void testDispVel(void);

    /// Test adding displacement and fault Lagrange multiplier subfields.
    void testDispLagrangeFault(void);

    /// Test adding pressure and fluid pressure subfields.
    void testPressure(void);

    /// Test adding displacement and temperature subfields.
    void testDispTemp(void);

    /// Test setValues().
    void testSetValues(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Initialze mesh, coordinate system, solution, and factory.
    void _initialize(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    SolutionFactory* _factory; ///< Test subject.
    TestSolutionFactory_Data* _data; ///< Test data.

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.
    pylith::topology::Field* _solution; ///< Solution field for test subject.

}; // class TestSolutionFactory

// ------------------------------------------------------------------------------------------------
class pylith::problems::TestSolutionFactory_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestSolutionFactory_Data(void);

    /// Destructor
    ~TestSolutionFactory_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    size_t dimension; ///< Spatial dimension.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    pylith::scales::Scales* scales; ///< Scales for nondimensionalization.

    std::map<std::string, pylith::topology::Field::SubfieldInfo> subfields;
    spatialdata::spatialdb::UserFunctionDB* solutionDB; ///< Spatial database with values for solution.

}; // class TestSolutionFactory_Data

// End of file
