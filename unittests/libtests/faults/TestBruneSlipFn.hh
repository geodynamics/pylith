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
 * @file unittests/libtests/faults/TestBruneSlipFn.hh
 *
 * @brief C++ TestBruneSlipFn object
 *
 * C++ unit testing for BruneSlipFn.
 */

#if !defined(pylith_faults_testbruneslipfn_hh)
#define pylith_faults_testbruneslipfn_hh

#include "pylith/utils/sievetypes.hh" // USES Mesh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestBruneSlipFn;
    class BruneSlipFn;

    namespace _TestBruneSlipFn {
      struct DataStruct;
    } // _BruneSlipTimeFn
  } // faults
} // pylith

/// C++ unit testing for BruneSlipFn
class pylith::faults::TestBruneSlipFn : public CppUnit::TestFixture
{ // class TestBruneSlipFn

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestBruneSlipFn );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testDbFinalSlip );
  CPPUNIT_TEST( testDbSlipTime );
  CPPUNIT_TEST( testDbPeakRate );
  CPPUNIT_TEST( testInitialize1D );
  CPPUNIT_TEST( testInitialize2D );
  CPPUNIT_TEST( testInitialize3D );
  CPPUNIT_TEST( testSlip );
  CPPUNIT_TEST( testSlipIncr );
  CPPUNIT_TEST( testSlipTH );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test dbFinalSlip().
  void testDbFinalSlip(void);

  /// Test dbSlipTime().
  void testDbSlipTime(void);

  /// Test dbPeakRate().
  void testDbPeakRate(void);

  /// Test initialize() in 1-D.
  void testInitialize1D(void);

  /// Test initialize() in 2-D.
  void testInitialize2D(void);

  /// Test initialize() in 3-D.
  void testInitialize3D(void);

  /// Test slip().
  void testSlip(void);

  /// Test slipIncr().
  void testSlipIncr(void);

  /// Test _slip().
  void testSlipTH(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize BruneSlipFn.
   *
   * @param faultMesh Fault mesh.
   * @param slipfn Brune slip function.
   */
  static
  void _initialize(ALE::Obj<Mesh>* faultMesh,
		   BruneSlipFn* slipfn);

  /** Test intialize().
   *
   * @param data Data for initialization and testing of BruneSlipFn.
   */
  static
  void _testInitialize(const _TestBruneSlipFn::DataStruct& data);

}; // class TestBruneSlipFn

#endif // pylith_faults_testbruneslipfn_hh


// End of file 
