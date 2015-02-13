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
 * @file unittests/libtests/faults/TestStepSlipFn.hh
 *
 * @brief C++ TestStepSlipFn object
 *
 * C++ unit testing for StepSlipFn.
 */

#if !defined(pylith_faults_teststepslipfn_hh)
#define pylith_faults_teststepslipfn_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/faults/faultsfwd.hh" // USES StepSlipFn
#include "pylith/topology/topologyfwd.hh" // USES Mesh

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestStepSlipFn;

    namespace _TestStepSlipFn {
      struct DataStruct;
    } // _BruneSlipTimeFn
  } // faults
} // pylith

/// C++ unit testing for StepSlipFn
class pylith::faults::TestStepSlipFn : public CppUnit::TestFixture
{ // class TestStepSlipFn

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestStepSlipFn );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testDbFinalSlip );
  CPPUNIT_TEST( testDbSlipTime );
  CPPUNIT_TEST( testInitialize2D );
  CPPUNIT_TEST( testInitialize3D );
  CPPUNIT_TEST( testSlip );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test dbFinalSlip().
  void testDbFinalSlip(void);

  /// Test dbSlipTime().
  void testDbSlipTime(void);

  /// Test initialize() in 2-D.
  void testInitialize2D(void);

  /// Test initialize() in 3-D.
  void testInitialize3D(void);

  /// Test slip().
  void testSlip(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize StepSlipFn.
   *
   * @param mesh Finite-element mesh of domain.
   * @param faultMesh Finite-element mesh of fault.
   * @param slipfn Step slip function.
   * @param originTime Origin time for earthquake rupture.
   */
  static
  void _initialize(topology::Mesh* mesh,
		   topology::Mesh* faultMesh,
		   StepSlipFn* slipfn,
		   const PylithScalar originTime);

  /** Test intialize().
   *
   * @param data Data for initialization and testing of StepSlipFn.
   */
  static
  void _testInitialize(const _TestStepSlipFn::DataStruct& data);

}; // class TestStepSlipFn

#endif // pylith_faults_teststepslipfn_hh


// End of file 
