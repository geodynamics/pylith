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
 * @file tests/libtests/materials/TestIsotropicLinearElasticity3D.hh
 *
 * @brief C++ TestIsotropicLinearElasticity3D object
 *
 * C++ unit testing for IsotropicLinearElasticity3D.
 */

#if !defined(pylith_materials_testisotropiclinearelasticity3d_hh)
#define pylith_materials_testisotropiclinearelasticity3d_hh

#include "TestMaterial.hh" // ISA TestMaterial

/// Namespace for pylith package
namespace pylith {
    namespace materials {
        class TestIsotropicLinearElasticity3D;

        class TestIsotropicLinearElasticity3D_Data;
    } // materials
} // pylith

/// C++ unit testing for IsotropicLinearElasticity3D
class pylith::materials::TestIsotropicLinearElasticity3D : public TestMaterial {

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticity3D, TestMaterial);

    // Tests specific to this materials parameters.
    CPPUNIT_TEST(testAccessors);

    // Tests that explicitly depend on details of this material.
    CPPUNIT_TEST(test_auxiliaryFieldSetup);
    CPPUNIT_TEST(testGetAuxField);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test useInertia(), useBodyForce(), useReferenceState().
    void testAccessors(void);

    /// Test _auxiliaryFieldSetup().
    void test_auxiliaryFieldSetup(void);

    /// Test getAuxField().
    void testGetAuxField(void);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get material.
     *
     * @returns Pointer to material.
     */
    Material* _material(void);

    /** Get test data.
     *
     * @returns Pointer to test data.
     */
    TestMaterial_Data* _data(void);

    /// Setup and populate solution fields.
    void _setupSolutionFields(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    IsotropicLinearElasticity3D* _mymaterial; ///< Object for testing.
    TestIsotropicLinearElasticity3D_Data* _mydata; ///< Data for testing.

}; // class TestIsotropicLinearElasticity3D

// =============================================================================
class pylith::materials::TestIsotropicLinearElasticity3D_Data : public TestMaterial_Data {

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    TestIsotropicLinearElasticity3D_Data(void);

    /// Destructor
    ~TestIsotropicLinearElasticity3D_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    // SPECIFIC TO MATERIAL, VALUES DEPEND ON TEST CASE
    double gravityVector[3]; ///< Array for gravity vector.

};

#endif // pylith_materials_testisotropiclinearelasticity3d_hh


// End of file
