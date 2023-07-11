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
 * @file tests/libtests/topology/TestFieldQuery.hh
 *
 * @brief C++ unit testing for FieldQuery.
 */

#if !defined(pylith_topology_testfieldquery_hh)
#define pylith_topology_testfieldquery_hh

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase::VectorFieldType

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA UserFunctionDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA CoordSys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional

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

    /// @defgroup Domain mesh information.
    /// @{
    int cellDim; ///< Cell dimension (matches space dimension).
    int numVertices; ///< Number of vertices.
    int numCells; ///< Number of cells.
    int numCorners; ///< Number of vertices per cell.
    const int* cells; ///< Array of vertices in cells [numCells*numCorners].
    const PylithScalar* coordinates; ///< Coordinates of vertices [numVertices*cellDim].

    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    spatialdata::units::Nondimensional* normalizer; ///< Scales for nondimensionalization.
    /// @}

    /// @defgroup Subfield discretization information
    /// @{
    int numAuxSubfields; ///< Number of auxiliary subfields.
    const char** auxSubfields; ///< Names of auxiliary subfields.
    pylith::topology::FieldBase::Description* auxDescriptions; ///< Descriptions for auxiliary subfields.
    pylith::topology::FieldBase::Discretization* auxDiscretizations; ///< Discretizations for auxiliary subfields.
    spatialdata::spatialdb::UserFunctionDB* auxDB; ///< Spatial database with auxiliary field.
    /// @}

}; // TestFieldQuery_Data

#endif // pylith_topology_testfieldquery_hh

// End of file
