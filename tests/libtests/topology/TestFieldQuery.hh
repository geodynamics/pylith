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

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase::VectorFieldType
#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA UserFunctionDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA CoordSys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Scales

namespace pylith {
    namespace topology {
        class TestFieldQuery;
        class TestFieldQuery_Data;
    } // topology
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestFieldQuery : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestFieldQuery(TestFieldQuery_Data* data);

    /// Destructor.
    ~TestFieldQuery(void);

    /// Test constructor.
    void testConstructor(void);

    /// Test setQuery().
    void testSetQuery(void);

    /// Test openDB(), closeDB().
    void testOpenClose(void);

    /// Test queryDB().
    void testQuery(void);

    /// Test queryDB() with NULL database.
    void testQueryNull(void);

    /// Test validatorPositive().
    void testValidatorPositive(void);

    /// Test validatorNonnegative().
    void testValidatorNonnegative(void);

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    /// Initialize mesh and test field.
    void _initialize(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    TestFieldQuery_Data* _data; ///< Data for testing.
    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.
    pylith::topology::Field* _field; ///< Field associated with mesh.
    pylith::topology::FieldQuery* _query; ///< Test field query associated with field.

    static const double FILL_VALUE; ///< Fill value for auxiliary field.

}; // class TestFieldQuery

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestFieldQuery_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestFieldQuery_Data(void);

    /// Destructor
    ~TestFieldQuery_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    pylith::meshio::MeshBuilder::Topology* topology; ///< Topology for domain mesh
    pylith::meshio::MeshBuilder::Geometry* geometry; ///< Geometry for domain mesh
    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    spatialdata::units::Scales* scales; ///< Scales for nondimensionalization.

    /// @defgroup Subfield discretization information
    /// @{
    int numAuxSubfields; ///< Number of auxiliary subfields.
    const char** auxSubfields; ///< Names of auxiliary subfields.
    pylith::topology::FieldBase::Description* auxDescriptions; ///< Descriptions for auxiliary subfields.
    pylith::topology::FieldBase::Discretization* auxDiscretizations; ///< Discretizations for auxiliary subfields.
    spatialdata::spatialdb::UserFunctionDB* auxDB; ///< Spatial database with auxiliary field.
    /// @}

}; // TestFieldQuery_Data

// End of file
