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

#if !defined(pylith_topology_testfieldmesh_hh)
#define pylith_topology_testfieldmesh_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/utils/petscfwd.h" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase::VectorFieldType

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
    namespace topology {
        class TestFieldMesh;
        class TestFieldMesh_Data;
    } // topology
} // pylith

// TestFieldMesh -------------------------------------------------------------
/// C++ unit testing for Field.
class pylith::topology::TestFieldMesh : public CppUnit::TestFixture {

    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE( TestFieldMesh );

    CPPUNIT_TEST( testConstructor );
    CPPUNIT_TEST( testMesh );
    CPPUNIT_TEST( testGeneralAccessors );
    CPPUNIT_TEST( testSectionAccessors );
    CPPUNIT_TEST( testVectorAccessors );
    CPPUNIT_TEST( testNewSection );
    CPPUNIT_TEST( testCloneSection );
    CPPUNIT_TEST( testSubfieldAccessors );
    CPPUNIT_TEST( testClear );
    CPPUNIT_TEST( testAllocate );
    CPPUNIT_TEST( testZeroLocal );
    CPPUNIT_TEST( testComplete );
    CPPUNIT_TEST( testCopy );
    CPPUNIT_TEST( testCopySubfield );
    CPPUNIT_TEST( testDimensionalize );
    CPPUNIT_TEST( testView );
    CPPUNIT_TEST( testScatter );

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS //////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test constructor.
    void testConstructor(void);

    /// Test mesh().
    void testMesh(void);

    /// Test label(), vectorFieldType(), scale(), addDimensionOkay(), spaceDim().
    void testGeneralAccessors(void);

    /// Test chartSize(), sectionSize(), localSection(), globalSection().
    void testSectionAccessors(void);

    /// Test localVector(), globalVector().
    void testVectorAccessors(void);

    /// Test newSection(points), newSection(domain), newSection(field).
    void testNewSection(void);

    /// Test cloneSection().
    void testCloneSection(void);

    /// Test subfieldAdd(), subfieldsSetup(), subfieldSetDof(), hasSubfield(), subfieldNames(), subfieldInfo().
    void testSubfieldAccessors(void);

    /// Test clear().
    void testClear(void);

    /// Test allocate().
    void testAllocate(void);

    /// Test zeroLocal().
    void testZeroLocal(void);

    /// Test complete().
    void testComplete(void);

    /// Test copy().
    void testCopy(void);

    /// Test copySubfield().
    void testCopySubfield(void);

    /// Test dimensionalize().
    void testDimensionalize(void);

    /// Test view().
    void testView(void);

    /// Test createScatter(), createScatterWithBC(), scatterLocalToGlobal(), scatterGlobalToLocal(), scatterVector().
    void testScatter(void);

    // PROTECTED METHODS ///////////////////////////////////////////////////////
protected:

    /// Initialize mesh and test field.
    void _initialize(void);

    /** Verify values in field match expected values.
     *
     * @param field Field containing values to test.
     * @param scale Scale to apply to expected values.
     */
    void _checkValues(const Field& field,
                      const PylithReal scale =1.0);

    /** Verify values in PETSc vector match expected values.
     *
     * @param vec PETSc vec containing values to test.
     * @param scale Scale to apply to expected values.
     */
    void _checkValues(const PetscVec& vec,
                      const PylithReal scale =1.0);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////
protected:

    TestFieldMesh_Data* _data; ///< Data for testing.
    Mesh* _mesh; ///< Finite-element mesh.
    Field* _field; ///< Test field associated with mesh.

}; // class TestFieldMesh

// TestFieldMesh_Data-----------------------------------------------------------
class pylith::topology::TestFieldMesh_Data {

    // PUBLIC METHODS //////////////////////////////////////////////////////////
public:

    /// Constructor
    TestFieldMesh_Data(void);

    /// Destructor
    ~TestFieldMesh_Data(void);

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

};  // TestFieldMesh_Data


#endif // pylith_topology_testfieldmesh_hh


// End of file
