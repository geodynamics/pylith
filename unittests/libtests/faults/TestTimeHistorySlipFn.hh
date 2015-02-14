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
 * @file unittests/libtests/faults/TestTimeHistorySlipFn.hh
 *
 * @brief C++ TestTimeHistorySlipFn object
 *
 * C++ unit testing for TimeHistorySlipFn.
 */

#if !defined(pylith_faults_testtimehistoryslipfn_hh)
#define pylith_faults_testtimehistoryslipfn_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/faults/faultsfwd.hh" // USES TimeHistorySlipFn
#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES TimeHistory

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestTimeHistorySlipFn;

    namespace _TestTimeHistorySlipFn {
      struct DataStruct;
    } // _TimeHistorySlipTimeFn
  } // faults
} // pylith

/// C++ unit testing for TimeHistorySlipFn
class pylith::faults::TestTimeHistorySlipFn : public CppUnit::TestFixture
{ // class TestTimeHistorySlipFn

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestTimeHistorySlipFn );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testDbAmplitude );
  CPPUNIT_TEST( testDbSlipTime );
  CPPUNIT_TEST( testDbTimeHistory );
  CPPUNIT_TEST( testInitialize2D );
  CPPUNIT_TEST( testInitialize3D );
  CPPUNIT_TEST( testSlip );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test dbAmplitude().
  void testDbAmplitude(void);

  /// Test dbSlipTime().
  void testDbSlipTime(void);

  /// Test dbTimeHistory().
  void testDbTimeHistory(void);

  /// Test initialize() in 2-D.
  void testInitialize2D(void);

  /// Test initialize() in 3-D.
  void testInitialize3D(void);

  /// Test slip().
  void testSlip(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize TimeHistorySlipFn.
   *
   * @param mesh Finite-element mesh of domain.
   * @param faultMesh Finite-element mesh of fault.
   * @param slipfn Step slip function.
   * @param th Time history database.
   * @param originTime Origin time for earthquake rupture.
   */
  static
  void _initialize(topology::Mesh* mesh,
		   topology::Mesh* faultMesh,
		   TimeHistorySlipFn* slipfn,
		   spatialdata::spatialdb::TimeHistory* th,
		   const PylithScalar originTime);

  /** Test intialize().
   *
   * @param data Data for initialization and testing of TimeHistorySlipFn.
   */
  static
  void _testInitialize(const _TestTimeHistorySlipFn::DataStruct& data);

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

}; // class TestTimeHistorySlipFn

#endif // pylith_faults_testtimehistoryslipfn_hh


// End of file 
