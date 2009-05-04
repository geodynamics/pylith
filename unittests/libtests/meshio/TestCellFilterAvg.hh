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
 * @file unittests/libtests/meshio/TestCellFilterAvg.hh
 *
 * @brief C++ TestCellFilterAvg object
 *
 * C++ unit testing for CellFilterAvg.
 */

#if !defined(pylith_meshio_testcellfilteravg_hh)
#define pylith_meshio_testcellfilteravg_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace meshio {
    class TestCellFilterAvg;
  } // meshio
} // pylith

// TestCellFilterAvg ----------------------------------------------------
class pylith::meshio::TestCellFilterAvg : public CppUnit::TestFixture
{ // class TestCellFilterAvg

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestCellFilterAvg );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testFilter );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test filter()
  void testFilter(void);

}; // class TestCellFilterAvg

#endif // pylith_meshio_testcellfilteravg_hh


// End of file 
