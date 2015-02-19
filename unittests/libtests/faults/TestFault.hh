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
  CPPUNIT_TEST( testEdge );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test id()
  void testID(void);

  /// Test label()
  void testLabel(void);

  /// Test edge()
  void testEdge(void);

}; // class TestFault

#endif // pylith_faults_testfault_hh

// End of file 
