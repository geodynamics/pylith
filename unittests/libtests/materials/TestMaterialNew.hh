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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/materials/TestMaterialNew.hh
 *
 * @brief C++ abstract base class for testing material objects.
 */

#if !defined(pylith_materials_testmaterialnew_hh)
#define pylith_materials_testmaterialnew_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/materials/materialsfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/topology/Field.hh" // HASA FieldBase::Discretization
#include "petscds.h" // USES PetscPointFunc, PetsPointJac
#include "spatialdata/spatialdb/SimpleGridDB.hh" // HOLDSA SimpleGridDB

/// Namespace for pylith package
namespace pylith {
    namespace materials {
        class TestMaterialNew;

        class TestMaterialNew_Data; // test data
    }   // materials
} // pylith

/// C++ abstract base class for testing material objects.
class pylith::materials::TestMaterialNew : public CppUnit::TestFixture { // class TestMaterialNew

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE( TestMaterialNew );

    CPPUNIT_TEST( testHasAuxField );
    CPPUNIT_TEST( testAuxFieldsDiscretization );
    CPPUNIT_TEST( testAuxFieldsDB );
    CPPUNIT_TEST( testNormalizer );

    CPPUNIT_TEST( testVerifyConfiguration );

    CPPUNIT_TEST( testDimension );
    CPPUNIT_TEST( testId );
    CPPUNIT_TEST( testLabel );
    CPPUNIT_TEST( testInitialize );

    CPPUNIT_TEST( testComputeResidual );
    CPPUNIT_TEST( testComputeRHSJacobian );
    CPPUNIT_TEST( testComputeLHSJacobianImplicit );
    CPPUNIT_TEST( testComputeLHSJacobianInverseExplicit );
    CPPUNIT_TEST( testUpdateStateVars );

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    virtual
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test hasAuxField().
    void testHasAuxField(void);

    /// Test auxFieldsDiscretization().
    void testAuxFieldsDiscretization(void);

    /// Test auxFieldsDB().
    void testAuxFieldsDB(void);

    /// Test normalizer().
    void testNormalizer(void);

    /// Test verifyConfiguration().
    void testVerifyConfiguration(void);

    /// Test checkConstraints().
    void testCheckConstraints(void);

    /// Test dimension().
    void testDimension(void);

    /// Test id().
    void testId(void);

    /// Test label().
    void testLabel(void);

    /// Test initialize().
    void testInitialize(void);

    /// Test computeRHSResidual(), computeLHSResidual().
    void testComputeResidual(void);

    /// Test computeRHSJacobian().
    void testComputeRHSJacobian(void);

    /// Test computeLHSJacobianImplicit().
    void testComputeLHSJacobianImplicit(void);

    /// Test computeLHSJacobianInverseExplicit().
    void testComputeLHSJacobianInverseExplicit(void);

    /// Test updateStateVars().
    void testUpdateStateVars(void);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get material.
     *
     * @returns Pointer to material.
     */
    virtual
    MaterialNew* _material(void);

    /** Get test data.
     *
     * @returns Pointer to test data.
     */
    virtual
    TestMaterialNew_Data* _data(void);

    /// Do minimal initilaization of test data.
    void _initializeMin(void);

    /// Do full initilaization of test data.
    void _initializeFull(void);

    /** Set field to zero on the boundary.
     *
     * @param[out] field Field in which to set boundary values to zero.
     */
    void _zeroBoundary(pylith::topology::Field* field);

    /** Setup and populate solution field.
     *
     * @param[out] field Solution field to setup and populate.
     * @param[in] dbFilename Filename for spatial database with values for field.
     * @param[in] isClone True if field is a clone (don't need full setup).
     */
    virtual
    void _setupSolutionField(pylith::topology::Field* field,
                             const char* dbFilename,
                             const bool isClone =false);


    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    // TestMaterialNew
    topology::Mesh* _mesh;   ///< Finite-element mesh.
    topology::Field* _solution1;   ///< Solution field.
    topology::Field* _solution2;   ///< Solution field.
    topology::Field* _solution1Dot;   ///< Time derivative of solution field.
    topology::Field* _solution2Dot;   ///< Time derivative of solution field.
    spatialdata::spatialdb::SimpleGridDB* _auxDB;   ///< Spatial database with data for auxiliary fields.

}; // class TestMaterialNew


// =============================================================================
class pylith::materials::TestMaterialNew_Data {

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    TestMaterialNew_Data(void);

    /// Destructor
    ~TestMaterialNew_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    int dimension; ///< Dimension of material.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    const char* materialLabel; ///< Label defining cells associated with material.
    int materialId; ///< Material id.
    const char* boundaryLabel; ///< Group defining domain boundary.

    PylithReal lengthScale; ///< Length scale for nondimensionalization.
    PylithReal timeScale; ///< Time scale for nondimensionalization.
    PylithReal pressureScale; ///< Pressure scale for nondimensionalization.
    PylithReal densityScale; ///< Density scale for nondimensionalization.

    PylithReal t; ///< Time for solution in simulation.
    PylithReal dt; ///< Time step in simulation.
    PylithReal tshift; ///< Time shift for LHS Jacobian.

    int numSolnFields; ///< Number of solution fields.
    pylith::topology::Field::DiscretizeInfo* solnDiscretizations; ///< Discretizations for solution fields.
    const char* solnDBFilename; ///< Name of file with data for solution.
    const char* pertDBFilename; ///< Name of file with data for perturbation.

    int numAuxFields; ///< Number of auxiliary fields.
    const char** auxFields; ///< Names of auxiliary fields.
    pylith::topology::Field::DiscretizeInfo* auxDiscretizations; ///< Discretizations for auxiliary fields.
    const char* auxDBFilename; ///< Name of file with data for auxFieldsDB.

    bool isExplicit; ///< True for explicit time stepping.
};


#endif // pylith_materials_testmaterialnew_hh


// End of file
