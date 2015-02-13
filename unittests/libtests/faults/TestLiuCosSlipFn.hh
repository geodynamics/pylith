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
 * @file unittests/libtests/faults/TestLiuCosSlipFn.hh
 *
 * @brief C++ TestLiuCosSlipFn object
 *
 * C++ unit testing for LiuCosSlipFn.
 */

#if !defined(pylith_faults_testliucosslipfn_hh)
#define pylith_faults_testliucosslipfn_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/faults/faultsfwd.hh" // USES LiuCosSlipFn
#include "pylith/topology/topologyfwd.hh" // USES Mesh

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestLiuCosSlipFn;

    namespace _TestLiuCosSlipFn {
      struct DataStruct;
    } // _LiuCosSlipTimeFn
  } // faults
} // pylith

/// C++ unit testing for LiuCosSlipFn
class pylith::faults::TestLiuCosSlipFn : public CppUnit::TestFixture
{ // class TestLiuCosSlipFn

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestLiuCosSlipFn );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testDbFinalSlip );
  CPPUNIT_TEST( testDbSlipTime );
  CPPUNIT_TEST( testDbRiseTime );
  CPPUNIT_TEST( testInitialize2D );
  CPPUNIT_TEST( testInitialize3D );
  CPPUNIT_TEST( testSlip );
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

  /// Test dbRiseTime().
  void testDbRiseTime(void);

  /// Test initialize() in 2-D.
  void testInitialize2D(void);

  /// Test initialize() in 3-D.
  void testInitialize3D(void);

  /// Test slip().
  void testSlip(void);

  /// Test _slip().
  void testSlipTH(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize LiuCosSlipFn.
   *
   * @param mesh Finite-element mesh of domain.
   * @param faultMesh Finite-element mesh of fault.
   * @param slipfn Step slip function.
   * @param originTime Origin time for earthquake rupture.
   */
  static
  void _initialize(topology::Mesh* mesh,
		   topology::Mesh* faultMesh,
		   LiuCosSlipFn* slipfn,
		   const PylithScalar originTime);

  /** Test intialize().
   *
   * @param data Data for initialization and testing of LiuCosSlipFn.
   */
  static
  void _testInitialize(const _TestLiuCosSlipFn::DataStruct& data);

  /** Slip time function.
   *
   * @param t Time relative to when slip begins.
   * @param finalSlip Final slip.
   * @param riseTime Rise time (t95).
   */
  static
  PylithScalar _slipFn(const PylithScalar t,
		 const PylithScalar finalSlip,
		 const PylithScalar riseTime);

}; // class TestLiuCosSlipFn

#endif // pylith_faults_testliucosslipfn_hh


// End of file 
