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

#include "pylith/faults/faultsfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // USES Mesh, SubMesh
#include "pylith/feassemble/feassemblefwd.hh" // HOLDSA Quadrature

#include <vector> // HASA std::vector

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKin;

    class CohesiveKinData;
  } // faults
} // pylith

/// C++ unit testing for FaultCohesiveKin
class pylith::faults::TestFaultCohesiveKin : public CppUnit::TestFixture
{ // class TestFaultCohesiveKin

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKin );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testEqsrc );
  CPPUNIT_TEST( testNeedNewJacobian );
  CPPUNIT_TEST( testUseSolnIncr );
  CPPUNIT_TEST( testUseLagrangeConstraints );

  CPPUNIT_TEST_SUITE_END();

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  CohesiveKinData* _data; ///< Data for testing
  feassemble::Quadrature<topology::SubMesh>* _quadrature; ///< Fault quad.
  std::vector<EqKinSrc*> _eqsrcs; ///< Array of Kinematic earthquake sources.
  std::vector<BruneSlipFn*> _slipfns; ///< Slip time function.
  bool _flipFault; ///< If true, flip fault orientation.

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

  /// Test needNewJacobian()
  void testNeedNewJacobian(void);

  /// Test useSolnIncr()
  void testUseSolnIncr(void);

  /// Test useLagrangeConstraints().
  void testUseLagrangeConstraints(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test integrateResidual().
  void testIntegrateResidual(void);

  /// Test integrateJacobian().
  void testIntegrateJacobian(void);

  /// Test integrateResidualAssembled().
  void testIntegrateResidualAssembled(void);

  /// Test integrateJacobianAssembled().
  void testIntegrateJacobianAssembled(void);

  /// Test updateState().
  void testUpdateState(void);

  /// Test _calcTractionsChange().
  void testCalcTractionsChange(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize FaultCohesiveKin interface condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param fault Cohesive fault interface condition to initialize.
   */
  void _initialize(topology::Mesh* mesh,
		   FaultCohesiveKin* const fault) const;

  /** Determine if vertex is a Lagrange multiplier constraint vertex.
   *
   * @param vertex Label of vertex.
   *
   * @returns True if vertex is a constraint vertex, false otherwise.
   */
  bool _isConstraintVertex(const int vertex) const;
  

}; // class TestFaultCohesiveKin

#endif // pylith_faults_testfaultcohesivekin_hh


// End of file 
