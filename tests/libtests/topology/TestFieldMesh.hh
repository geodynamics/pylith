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
 * @file tests/libtests/topology/TestFieldMesh.hh
 *
 * @brief C++ unit testing for Field.
 */

#if !defined(pylith_topology_testfieldmesh_hh)
#define pylith_topology_testfieldmesh_hh

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/utils/petscfwd.h" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Description

namespace pylith {
    namespace topology {
        class TestFieldMesh;
        class TestFieldMesh_Data;
    } // topology
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestFieldMesh : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestFieldMesh(TestFieldMesh_Data* data);

    /// Destructor.
    ~TestFieldMesh(void);

    /// Test constructor.
    void testConstructor(void);

    /// Test copy constructor.
    void testCopyConstructor(void);

    /// Test mesh().
    void testMesh(void);

    /// Test getLabel(), vectorFieldType(), scale(), addDimensionOkay(), getSpaceDim().
    void testGeneralAccessors(void);

    /// Test chartSize(), getStorageSize(), selectedSection(), globalSection().
    void testSectionAccessors(void);

    /// Test localVector(), globalVector().
    void testVectorAccessors(void);

    /// Test subfieldAdd(), subfieldsSetup(), hasSubfield(), subfieldNames(), subfieldInfo().
    void testSubfieldAccessors(void);

    /// Test allocate().
    void testAllocate(void);

    /// Test zeroLocal().
    void testZeroLocal(void);

    /// Test view().
    void testView(void);

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    /// Initialize mesh and test field.
    void _initialize(void);

    /** Verify values in field match expected values.
     *
     * @param field Field containing values to test.
     * @param scale Scale to apply to expected values.
     */
    void _checkValues(const Field& field,
                      const PylithReal scale=1.0);

    /** Verify values in PETSc vector match expected values.
     *
     * @param vec PETSc vec containing values to test.
     * @param scale Scale to apply to expected values.
     */
    void _checkValues(const PetscVec& vec,
                      const PylithReal scale=1.0);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////
protected:

    TestFieldMesh_Data* _data; ///< Data for testing.
    Mesh* _mesh; ///< Finite-element mesh.
    Field* _field; ///< Test field associated with mesh.

}; // class TestFieldMesh

// ------------------------------------------------------------------------------------------------
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
    int numCells; ///< Number of cells.
    int numCorners; ///< Number of vertices per cell.
    int* cells; ///< Array of vertices in cells [numCells*numCorners].
    PylithScalar* coordinates; ///< Coordinates of vertices [numVertices*cellDim].
    /// @}

    /// @defgroup Subfield A information.
    /// @{
    pylith::topology::FieldBase::Description descriptionA;
    pylith::topology::FieldBase::Discretization discretizationA;
    PylithScalar* subfieldAValues; ///< Array of values to for subfield,

    const char* bcALabel; ///< Label for boundary condition.
    int bcALabelId; ///< Label id for boundary condition.
    int bcANumConstrainedDOF; ///< Number of constrained DOF for boundary condition.
    int* bcAConstrainedDOF; ///< Array of constrained DOF.
    int bcANumVertices; ///< Number of vertices associated with boundary condition.
    int* bcAVertices; ///< Array of vertex indices.
    /// @}

    /// @defgroup Subfield B information.
    /// @{
    pylith::topology::FieldBase::Description descriptionB;
    pylith::topology::FieldBase::Discretization discretizationB;
    PylithScalar* subfieldBValues; ///< Array of values to for subfield,

    const char* bcBLabel; ///< Label for boundary condition.
    int bcBLabelId; ///< Label id for boundary condition.
    int bcBNumConstrainedDOF; ///< Number of constrained DOF for boundary condition.
    int* bcBConstrainedDOF; ///< Array of constrained DOF.
    int bcBNumVertices; ///< Number of vertices associated with boundary condition.
    int* bcBVertices; ///< Array of vertex indices.
    /// @}

}; // TestFieldMesh_Data

#endif // pylith_topology_testfieldmesh_hh

// End of file
