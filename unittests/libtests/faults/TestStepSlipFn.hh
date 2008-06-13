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
 * @file unittests/libtests/faults/TestStepSlipFn.hh
 *
 * @brief C++ TestStepSlipFn object
 *
 * C++ unit testing for StepSlipFn.
 */

#if !defined(pylith_faults_teststepslipfn_hh)
#define pylith_faults_teststepslipfn_hh

#include "pylith/utils/sievetypes.hh" // USES Mesh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestStepSlipFn;
    class StepSlipFn;

    namespace _TestStepSlipFn {
      struct DataStruct;
    } // _StepSlipTimeFn
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
  CPPUNIT_TEST( testInitialize1D );
  CPPUNIT_TEST( testInitialize2D );
  CPPUNIT_TEST( testInitialize3D );
  CPPUNIT_TEST( testSlip );
  CPPUNIT_TEST( testSlipIncr );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test dbFinalSlip().
  void testDbFinalSlip(void);

  /// Test dbSlipTime().
  void testDbSlipTime(void);

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

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize StepSlipFn.
   *
   * @param faultMesh Fault mesh.
   * @param slipfn Step slip function.
   * @param originTime Origin time for earthquake rupture.
   */
  static
  void _initialize(ALE::Obj<Mesh>* faultMesh,
		   StepSlipFn* slipfn,
		   const double originTime);

  /** Test intialize().
   *
   * @param data Data for initialization and testing of StepSlipFn.
   */
  static
  void _testInitialize(const _TestStepSlipFn::DataStruct& data);

}; // class TestStepSlipFn

#endif // pylith_faults_teststepslipfn_hh


// End of file 
