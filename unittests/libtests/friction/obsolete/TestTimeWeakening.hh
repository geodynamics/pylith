// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/friction/TestTimeWeakening.hh
 *
 * @brief C++ TestTimeWeakening object
 *
 * C++ unit testing for TimeWeakening.
 */

#if !defined(pylith_friction_testtimeweakening_hh)
#define pylith_friction_testtimeweakening_hh

#include "TestFrictionModel.hh"

/// Namespace for pylith package
namespace pylith {
  namespace friction {
    class TestTimeWeakening;
  } // friction
} // pylith

/// C++ unit testing for TimeWeakening
class pylith::friction::TestTimeWeakening : public TestFrictionModel
{ // class TestTimeWeakening

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestTimeWeakening );

  CPPUNIT_TEST( testPropertiesMetadata );
  CPPUNIT_TEST( testStateVarsMetadata );
  CPPUNIT_TEST( testDBToProperties );
  CPPUNIT_TEST( testNonDimProperties );
  CPPUNIT_TEST( testDimProperties );
  CPPUNIT_TEST( testDBToStateVars );
  CPPUNIT_TEST( testNonDimStateVars );
  CPPUNIT_TEST( testDimStateVars );
  CPPUNIT_TEST( testHasPropStateVar );
  CPPUNIT_TEST( test_calcFriction );
  CPPUNIT_TEST( test_calcFrictionDeriv );
  CPPUNIT_TEST( test_updateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Test properties metadata.
  void testPropertiesMetadata(void);

  /// Test state variable metadata.
  void testStateVarsMetadata(void);

  /// Test hasPropStateVar().
  void testHasPropStateVar(void);

}; // class TestTimeWeakening

#endif // pylith_friction_testtimeweakening_hh


// End of file 
