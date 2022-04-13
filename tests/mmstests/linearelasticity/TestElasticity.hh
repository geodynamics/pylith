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

/**
 * @file mmstets/elasticity/TestElasticity.hh
 *
 * @brief C++ abstract base class for testing Elasticity for various rheologies.
 */

#if !defined(pylith_mmstests_testelasticity_hh)
#define pylith_mmstests_testelasticity_hh

#include "pylith/testing/MMSTest.hh" // ISA MMSTEST

#include "pylith/materials/materialsfwd.hh" // forward declarations
#include "pylith/bc/bcfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/topology/Field.hh" // HASA FieldBase::Discretization

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA UserFunctionDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA CoordSys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional

namespace pylith {
    namespace mmstests {
        class TestElasticity;

        class TestElasticity_Data; // test data
    } // mmstets
} // pylith

/// C++ abstract base class for testing Elasticity with various rheologies.
class pylith::mmstests::TestElasticity : public pylith::testing::MMSTest {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    virtual
    void setUp(void);

    /// Deallocate testing data.
    virtual
    void tearDown(void);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Initialize objects for test.
    void _initialize(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::materials::Elasticity* _material; ///< Material for testing.
    pylith::bc::DirichletUserFn* _bc; ///< Dirichlet boundary condition.
    TestElasticity_Data* _data; ///< Test parameters.

}; // class TestElasticity

// =====================================================================================================================
class pylith::mmstests::TestElasticity_Data {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestElasticity_Data(void);

    /// Destructor
    ~TestElasticity_Data(void);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    int spaceDim; ///< Spatial dimension of domain.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    const char* boundaryLabel; ///< Group defining domain boundary.

    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    spatialdata::spatialdb::GravityField* gravityField; ///< Gravity field.
    spatialdata::units::Nondimensional* normalizer; ///< Scales for nondimensionalization.

    PylithReal startTime; ///< Start time for simulation.
    PylithReal endTime; ///< End time for simulation.
    PylithReal timeStep; ///< Time step in simulation.

    int numSolnSubfields; ///< Number of solution fields.
    pylith::topology::Field::Discretization* solnDiscretizations; ///< Discretizations for solution fields.

    int numAuxSubfields; ///< Number of auxiliary subfields.
    const char** auxSubfields; ///< Names of auxiliary subfields.
    pylith::topology::Field::Discretization* auxDiscretizations; ///< Discretizations for auxiliary subfields.
    spatialdata::spatialdb::UserFunctionDB* auxDB; ///< Spatial database with auxiliary field.

    bool isExplicit; ///< True for explicit time stepping.
}; // TestElasticity_Data

#endif // pylith_mmstests_testelasticity_hh

// End of file
