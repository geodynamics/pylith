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
 * @file unittests/libtests/friction/TestRateStateAgeing.hh
 *
 * @brief C++ TestRateStateAgeing object
 *
 * C++ unit testing for RateStateAgeing.
 */

#if !defined(pylith_friction_testelasticisotropic3d_hh)
#define pylith_friction_testelasticisotropic3d_hh

#include "TestFrictionModel.hh"

/// Namespace for pylith package
namespace pylith {
  namespace friction {
    class TestRateStateAgeing;
  } // friction
} // pylith

/// C++ unit testing for RateStateAgeing
class pylith::friction::TestRateStateAgeing : public TestFrictionModel
{ // class TestRateStateAgeing

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestRateStateAgeing );

  CPPUNIT_TEST( testPropertiesMetadata );
  CPPUNIT_TEST( testStateVarsMetadata );
  CPPUNIT_TEST( testDBToProperties );
  CPPUNIT_TEST( testNonDimProperties );
  CPPUNIT_TEST( testDimProperties );
  CPPUNIT_TEST( testDBToStateVars );
  CPPUNIT_TEST( testNonDimStateVars );
  CPPUNIT_TEST( testDimStateVars );
  CPPUNIT_TEST( testHasProperty );
  CPPUNIT_TEST( testHasStateVar );
  CPPUNIT_TEST( test_calcFriction );
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

  /// Test hasProperty().
  void testHasProperty(void);

  /// Test hasStateVar().
  void testHasStateVar(void);

}; // class TestRateStateAgeing

#endif // pylith_friction_testelasticisotropic3d_hh


// End of file 
