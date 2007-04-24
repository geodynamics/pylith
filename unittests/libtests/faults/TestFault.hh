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
 * @file unittests/libtests/faults/TestFault.hh
 *
 * @brief C++ TestFault object
 *
 * C++ unit testing for Fault.
 */

#if !defined(pylith_faults_testfault_hh)
#define pylith_faults_testfault_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class Fault;
    class TestFault;
  } // faults
} // pylith

/// C++ unit testing for Fault
class pylith::faults::TestFault : public CppUnit::TestFixture
{ // class TestFault

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFault );
  CPPUNIT_TEST( testID );
  CPPUNIT_TEST( testLabel );
#if 0
  CPPUNIT_TEST( testInitialize );
#endif
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test id()
  void testID(void);

  /// Test label()
  void testLabel(void);

#if 0
  /// Test initialize()
  void testInitialize(void);
#endif

}; // class TestFault

#endif // pylith_faults_testfault_hh

// End of file 
