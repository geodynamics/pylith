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
    class FaultCohesiveKinData; // HOLDSA FaultCohesiveKinData
  } // faults
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

  FaultCohesiveKinData* _data; ///< Data for testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

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

}; // class TestFaultCohesiveKin

#endif // pylith_faults_testfaultcohesivekin_hh


// End of file 
