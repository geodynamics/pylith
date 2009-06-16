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
 * @file unittests/libtests/bc/TestTimeDependent.hh
 *
 * @brief C++ TestTimeDependent object.
 *
 * C++ unit testing for TimeDependent.
 */

#if !defined(pylith_bc_testtimedependent_hh)
#define pylith_bc_testtimedependent_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestTimeDependent;
  } // bc
} // pylith

/// C++ unit testing for PointForce.
class pylith::bc::TestTimeDependent : public CppUnit::TestFixture
{ // class TestTimeDependent

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestTimeDependent );

  CPPUNIT_TEST( testDBInitial );
  CPPUNIT_TEST( testDBRate );
  CPPUNIT_TEST( testDBChange );
  CPPUNIT_TEST( testDBTimeHistory );
  CPPUNIT_TEST( testVerifyConfiguration );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test dbInitial().
  void testDBInitial(void);

  /// Test dbRate().
  void testDBRate(void);

  /// Test dbChange().
  void testDBChange(void);

  /// Test dbTimeHistory().
  void testDBTimeHistory(void);

  /// Test verifyConfiguration().
  void testVerifyConfiguration(void);

}; // class TestTimeDependent

#endif // pylith_bc_pointforce_hh


// End of file 
