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
 * @file unittests/libtests/faults/TestEqKinSrc.hh
 *
 * @brief C++ TestEqKinSrc object
 *
 * C++ unit testing for EqKinSrc.
 */

#if !defined(pylith_faults_testeqkinsrc_hh)
#define pylith_faults_testeqkinsrc_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestEqKinSrc;

    namespace _TestEqKinSrc {
      struct DataStruct;
    } // _BruneSlipTimeFn
  } // faults
} // pylith

/// C++ unit testing for EqKinSrc
class pylith::faults::TestEqKinSrc : public CppUnit::TestFixture
{ // class TestEqKinSrc

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestEqKinSrc );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testSlipFn );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testSlip );
  CPPUNIT_TEST( testSlipIncr );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test slipFn().
  void testSlipFn(void);

  /// Test initialize(). Use 2-D mesh with Brune slip function to test
  /// initialize().
  void testInitialize(void);

  /// Test slip(). Use 2-D mesh with Brune slip function to test
  /// slip().
  void testSlip(void);

  /// Test slipIncr(). Use 2-D mesh with Brune slip function to test
  /// slip().
  void testSlipIncr(void);

}; // class TestEqKinSrc

#endif // pylith_faults_testeqkinsrc_hh


// End of file 
