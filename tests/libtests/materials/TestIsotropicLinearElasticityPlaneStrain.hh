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
 * @file tests/libtests/materials/TestIsotropicLinearElasticityPlaneStrain.hh
 *
 * @brief C++ TestIsotropicLinearElasticityPlaneStrain object
 *
 * C++ unit testing for IsotropicLinearElasticityPlaneStrain.
 */

#if !defined(pylith_materials_testisotropiclinearelasticityplanestrain_hh)
#define pylith_materials_testisotropiclinearelasticityplanestrain_hh

#include "TestMaterial.hh" // ISA TestMaterial

/// Namespace for pylith package
namespace pylith {
    namespace materials {
        class TestIsotropicLinearElasticityPlaneStrain;

        class TestIsotropicLinearElasticityPlaneStrain_Data;
    } // materials
} // pylith

/// C++ unit testing for IsotropicLinearElasticityPlaneStrain
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain : public TestMaterial {

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain, TestMaterial);

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

    IsotropicLinearElasticityPlaneStrain* _mymaterial; ///< Object for testing.
    TestIsotropicLinearElasticityPlaneStrain_Data* _mydata; ///< Data for testing.

}; // class TestIsotropicLinearElasticityPlaneStrain

// =============================================================================
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_Data : public TestMaterial_Data {

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    TestIsotropicLinearElasticityPlaneStrain_Data(void);

    /// Destructor
    ~TestIsotropicLinearElasticityPlaneStrain_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    // SPECIFIC TO MATERIAL, VALUES DEPEND ON TEST CASE
    double gravityVector[3]; ///< Array for gravity vector.

};

#endif // pylith_materials_testisotropiclinearelasticityplanestrain_hh


// End of file
