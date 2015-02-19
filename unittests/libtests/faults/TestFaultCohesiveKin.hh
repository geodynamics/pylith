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
#include "pylith/topology/topologyfwd.hh" // USES Mesh
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
  CPPUNIT_TEST( testUseLagrangeConstraints );

  CPPUNIT_TEST_SUITE_END();

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  CohesiveKinData* _data; ///< Data for testing
  feassemble::Quadrature* _quadrature; ///< Fault quadrature.a
  std::vector<EqKinSrc*> _eqsrcs; ///< Array of Kinematic earthquake sources.
  std::vector<BruneSlipFn*> _slipfns; ///< Slip time function.

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

  /// Test useLagrangeConstraints().
  void testUseLagrangeConstraints(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test integrateResidual().
  void testIntegrateResidual(void);

  /// Test integrateJacobian().
  void testIntegrateJacobian(void);

  /// Test integrateJacobian() with lumped Jacobian.
  void testIntegrateJacobianLumped(void);

  /// Test adjustSolnLumped().
  void testAdjustSolnLumped(void);

  /// Test _calcTractionsChange().
  void testCalcTractionsChange(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Set field values.
   *
   * @param field Field to assign values to.
   * @param data Array of data.
   * @param scale Scale to nondimensionalize values by.
   */
  void _fieldSetValues(topology::Field* field,
		       const PylithScalar* data,
		       const PylithScalar scale =1.0);

  /** Initialize FaultCohesiveKin interface condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param fault Cohesive fault interface condition to initialize.
   * @param fields Solution fields.
   */
  void _initialize(topology::Mesh* const mesh,
		   FaultCohesiveKin* const fault,
		   topology::SolutionFields* const fields) const;

  /** Determine if point is a Lagrange multiplier constraint point.
   *
   * @param point Label of point.
   *
   * @returns True if point is a constraint point, false otherwise.
   */
  bool _isConstraintEdge(const int point) const;
  

}; // class TestFaultCohesiveKin

#endif // pylith_faults_testfaultcohesivekin_hh


// End of file 
