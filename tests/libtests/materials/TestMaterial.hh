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
 * @file tests/libtests/materials/TestMaterial.hh
 *
 * @brief C++ abstract base class for testing material objects.
 */

#if !defined(pylith_materials_testMaterial_hh)
#define pylith_materials_testMaterial_hh

#include <cppunit/extensions/HelperMacros.h>
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/materials/materialsfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/topology/Field.hh" // HASA FieldBase::Discretization

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA UserFunctionDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA CoordSys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional

/// Namespace for pylith package
namespace pylith {
    namespace materials {
        class TestMaterial;

        class TestMaterial_Data; // test data
    } // materials
} // pylith

/// C++ abstract base class for testing material objects.
class pylith::materials::TestMaterial : public CppUnit::TestFixture, public pylith::utils::GenericComponent {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestMaterial);

    CPPUNIT_TEST(testAuxField);
    CPPUNIT_TEST(testAuxSubfieldDiscretization);
    CPPUNIT_TEST(testAuxFieldDB);
    CPPUNIT_TEST(testNormalizer);

    CPPUNIT_TEST(testVerifyConfiguration);

    CPPUNIT_TEST(testAccessors);
    CPPUNIT_TEST(testInitialize);

    CPPUNIT_TEST(testComputeResidual);
    CPPUNIT_TEST(testComputeJacobian);
    CPPUNIT_TEST(testComputeLHSJacobianImplicit);
    CPPUNIT_TEST(testComputeLHSJacobianInverseExplicit);
    CPPUNIT_TEST(testUpdateStateVars);

    CPPUNIT_TEST_SUITE_END_ABSTRACT();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    virtual
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test auxField().
    void testAuxField(void);

    /// Test auxSubfieldDiscretization().
    void testAuxSubfieldDiscretization(void);

    /// Test auxFieldDB().
    void testAuxFieldDB(void);

    /// Test normalizer().
    void testNormalizer(void);

    /// Test verifyConfiguration().
    void testVerifyConfiguration(void);

    /// Test checkConstraints().
    void testCheckConstraints(void);

    /// Test dimension(), id(), and getLabel().
    void testAccessors(void);

    /// Test initialize().
    void testInitialize(void);

    /// Test computeRHSResidual(), computeLHSResidual().
    void testComputeResidual(void);

    /// Test computeJacobian().
    void testComputeJacobian(void);

    /// Test computeLHSJacobianImplicit().
    void testComputeLHSJacobianImplicit(void);

    /// Test computeLHSJacobianInverseExplicit().
    void testComputeLHSJacobianInverseExplicit(void);

    /// Test _updateStateVars().
    void testUpdateStateVars(void);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get material.
     *
     * @returns Pointer to material.
     */
    virtual
    Material* _material(void) = 0;

    /** Get test data.
     *
     * @returns Pointer to test data.
     */
    virtual
    TestMaterial_Data* _data(void) = 0;

    /// Do minimal initilaization of test data.
    void _initializeMin(void);

    /// Do full initilaization of test data.
    void _initializeFull(void);

    /** Set field (and, optionally, matrix rows and columns) to zero on the boundary.
     *
     * @param[out] field Field in which to set boundary values to zero.
     */
    void _zeroBoundary(pylith::topology::Field* field);

    /// Setup and populate solution fields.
    virtual
    void _setupSolutionFields(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    // TestMaterial
    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.
    pylith::topology::Fields* _solutionFields; ///< Contrainer for solution fields.

}; // class TestMaterial

// =============================================================================
class pylith::materials::TestMaterial_Data {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    TestMaterial_Data(void);

    /// Destructor
    ~TestMaterial_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    int dimension; ///< Dimension of material.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    const char* boundaryLabel; ///< Group defining domain boundary.

    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    spatialdata::spatialdb::GravityField* gravityField; ///< Gravity field.
    spatialdata::units::Nondimensional* normalizer; ///< Scales for nondimensionalization.

    PylithReal t; ///< Time for solution in simulation.
    PylithReal dt; ///< Time step in simulation.
    PylithReal s_tshift; ///< Time shift for LHS Jacobian.
    PylithReal perturbation; ///< Maximum amplitude of random perturbation.

    int numSolnSubfields; ///< Number of solution fields.
    pylith::topology::Field::Discretization* solnDiscretizations; ///< Discretizations for solution fields.
    spatialdata::spatialdb::UserFunctionDB* solnDB; ///< Spatial database with solution.
    spatialdata::spatialdb::UserFunctionDB* perturbDB; ///< Spatial database with solution + perturbation.

    int numAuxSubfields; ///< Number of auxiliary subfields.
    const char** auxSubfields; ///< Names of auxiliary subfields.
    pylith::topology::Field::Discretization* auxDiscretizations; ///< Discretizations for auxiliary subfields.
    spatialdata::spatialdb::UserFunctionDB* auxDB; ///< Spatial database with auxiliary field.
    spatialdata::spatialdb::UserFunctionDB* auxUpdateDB; ///< Spatial database with updated auxiliary field.

    bool isExplicit; ///< True for explicit time stepping.
};

#endif // pylith_materials_testMaterial_hh

// End of file
