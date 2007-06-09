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
 * @file unittests/libtests/faults/TestFaultCohesiveKin.hh
 *
 * @brief C++ TestFaultCohesiveKin object
 *
 * C++ unit testing for FaultCohesiveKin.
 */

#if !defined(pylith_faults_testfaultcohesivekin_hh)
#define pylith_faults_testfaultcohesivekin_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKin;

    class FaultCohesiveKin; // USES FaultCohesiveKin
    class CohesiveKinData; // HOLDSA CohesiveKinData
    class EqKinSrc; // HOLDSA EqKinSrc
    class BruneSlipFn; // HOLDSA BruneSlipFn
  } // faults

  namespace feassemble {
    class Quadrature; // HOLDSA Quadrature
  } // feassemble
} // pylith

/// C++ unit testing for FaultCohesiveKin
class pylith::faults::TestFaultCohesiveKin : public CppUnit::TestFixture
{ // class TestFaultCohesiveKin

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKin );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testEqsrc );
  CPPUNIT_TEST( testUseLagrangeConstraints );

  CPPUNIT_TEST_SUITE_END();

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  CohesiveKinData* _data; ///< Data for testing
  feassemble::Quadrature* _quadrature; ///< Data used in testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Test constructor.
  void testConstructor(void);

  /// Test eqsrc().
  void testEqsrc(void);

  /// Test useLagrangeConstraints().
  void testUseLagrangeConstraints(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test integrateResidual().
  void testIntegrateResidual(void);

  /// Test integrateJacobian().
  void testIntegrateJacobian(void);

  /// Test setConstraintSizes().
  void testSetConstraintSizes(void);

  /// Test setField().
  void testSetField(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize FaultCohesiveKin interface condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param fault Cohesive fault interface condition to initialize.
   */
  void _initialize(ALE::Obj<ALE::Mesh>* mesh,
		   FaultCohesiveKin* const fault) const;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  EqKinSrc* _eqsrc; ///< Kinematic earthquake source
  BruneSlipFn* _slipfn; ///< Slip time function

}; // class TestFaultCohesiveKin

#endif // pylith_faults_testfaultcohesivekin_hh


// End of file 
