// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/topology/TestFieldMesh.hh
 *
 * @brief C++ unit testing for Field.
 */

#if !defined(pylith_topology_testfieldquery_hh)
#define pylith_topology_testfieldquery_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase::VectorFieldType

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
    namespace topology {
        class TestFieldQuery;
        class TestFieldQuery_Data;
    } // topology
} // pylith

// TestFieldQuery -------------------------------------------------------------
/// C++ unit testing for Field.
class pylith::topology::TestFieldQuery : public CppUnit::TestFixture {

    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE( TestFieldQuery );

    CPPUNIT_TEST( testConstructor );
    CPPUNIT_TEST( testQueryFn );
    CPPUNIT_TEST( testOpenClose );
    CPPUNIT_TEST( testQuery );
    CPPUNIT_TEST( testDBQueryGeneric );
    CPPUNIT_TEST( testValidatorPositive );

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS //////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test constructor.
    void testConstructor(void);

    /// Test queryFn().
    void testQueryFn(void);

    /// Test openDB(), closeDB()..
    void testOpenClose(void);

    /// Test queryDB().
    void testQuery(void);

    /// Test dbQueryGeneric.
    void testDBQueryGeneric(void);

    /// Test validatorPositive().
    void testValidatorPositive(void);

    // PRIVATE METHODS /////////////////////////////////////////////////////////
private:

    /// Initialize mesh and test field.
    void _initialize(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////
protected:

    TestFieldQuery_Data* _data; ///< Data for testing.
    Mesh* _mesh; ///< Finite-element mesh.
    Field* _field; ///< Field associated with mesh.
    FieldQuery* _query; ///< Test field query associated with field.

}; // class TestFieldQuery

// TestFieldQuery_Data-----------------------------------------------------------
class pylith::topology::TestFieldQuery_Data {

    // PUBLIC METHODS //////////////////////////////////////////////////////////
public:

    /// Constructor
    TestFieldQuery_Data(void);

    /// Destructor
    ~TestFieldQuery_Data(void);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////
public:

    /// @defgroup Domain mesh information.
    /// @{
    int cellDim; ///< Cell dimension (matches space dimension).
    int numVertices; ///< Number of vertices.
    int numCells;   ///< Number of cells.
    int numCorners; ///< Number of vertices per cell.
    int* cells; ///< Array of vertices in cells [numCells*numCorners].
    PylithScalar* coordinates;  ///< Coordinates of vertices [numVertices*cellDim].
    /// @}

    /// @defgroup Subfield A information.
    /// @{
    const char* subfieldAName; ///< Name of subfield.
    FieldBase::VectorFieldEnum subfieldAType; ///< Vector field type.
    PylithReal subfieldAScale; ///< Scale.
    int subfieldANumComponents; ///< Number of components in subfield.
    const char** subfieldAComponents; ///< Names of components.
    PylithScalar* subfieldAValues; ///< Array of values to for subfield,
    int subfieldABasisOrder; ///< Basis order for discretization.
    int subfieldAQuadOrder; ///< Quadrature order for discretization.

    const char* bcALabel; ///< Label for boundary condition.
    int bcALabelId; ///< Label id for boundary condition.
    int bcANumConstrainedDOF; ///< Number of constrained DOF for boundary condition.
    int* bcAConstrainedDOF; ///< Array of constrained DOF.
    int bcANumVertices; ///< Number of vertices assocaited with boundary condition.
    int* bcAVertices; ///< Array of vertex indices.
    /// @}

    /// @defgroup Subfield B information.
    /// @{
    const char* subfieldBName; ///< Name of subfield.
    FieldBase::VectorFieldEnum subfieldBType; ///< Vector field type.
    PylithReal subfieldBScale; ///< Scale.
    int subfieldBNumComponents; ///< Number of components in subfield.
    const char** subfieldBComponents; ///< Names of components.
    PylithScalar* subfieldBValues; ///< Array of values to for subfield,
    int subfieldBBasisOrder; ///< Basis order for discretization.
    int subfieldBQuadOrder; ///< Quadrature order for discretization.

    const char* bcBLabel; ///< Label for boundary condition.
    int bcBLabelId; ///< Label id for boundary condition.
    int bcBNumConstrainedDOF; ///< Number of constrained DOF for boundary condition.
    int* bcBConstrainedDOF; ///< Array of constrained DOF.
    int bcBNumVertices; ///< Number of vertices assocaited with boundary condition.
    int* bcBVertices; ///< Array of vertex indices.
    /// @}

};  // TestFieldQuery_Data


#endif // pylith_topology_testfieldquery_hh


// End of file
