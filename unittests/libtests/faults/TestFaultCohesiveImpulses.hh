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
 * @file unittests/libtests/faults/TestFaultCohesiveImpulses.hh
 *
 * @brief C++ TestFaultCohesiveImpulses object
 *
 * C++ unit testing for FaultCohesiveImpulses.
 */

#if !defined(pylith_faults_testfaultcohesiveimpulses_hh)
#define pylith_faults_testfaultcohesiveimpulses_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/faults/faultsfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "pylith/feassemble/feassemblefwd.hh" // HOLDSA Quadrature
#include "pylith/friction/frictionfwd.hh" // HOLDSA FrictionModel
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB
#include <vector> // HASA std::vector
/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveImpulses;

    class CohesiveImpulsesData;
  } // faults
} // pylith

/// C++ unit testing for FaultCohesiveImpulses
class pylith::faults::TestFaultCohesiveImpulses: public CppUnit::TestFixture
{ // class TestFaultCohesiveImpulses

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveImpulses );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testDBImpulseAmp );
  CPPUNIT_TEST( testThreshold );
  CPPUNIT_TEST( testImpulseDOF );

  // Tests in derived classes:
  // testNumImpulses()
  // testInitialize()
  // testIntegrateResidual()

  CPPUNIT_TEST_SUITE_END();

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

  CohesiveImpulsesData* _data; ///< Data for testing
  feassemble::Quadrature* _quadrature; ///< Fault quadrature.
  spatialdata::spatialdb::SpatialDB* _dbImpulseAmp; ///< Initial tractions.

  // PUBLIC METHODS /////////////////////////////////////////////////////
public:

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Test constructor.
  void testConstructor(void);

  /// Test dbInitialTract().
  void testDBImpulseAmp(void);

  /// Test threshold().
  void testThreshold(void);

  /// Test impulseDOF().
  void testImpulseDOF(void);

  /// Test numImpulses().
  void testNumImpulses(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test integrateResidual().
  void testIntegrateResidual(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private:

  /** Initialize FaultCohesiveImpulses interface condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param fault Cohesive fault interface condition to initialize.
   * @param fields Solution fields.
   */
  void _initialize(topology::Mesh* const mesh,
		   FaultCohesiveImpulses* const fault,
		   topology::SolutionFields* const fields);

  /** Determine if point is a Lagrange multiplier constraint point.
   *
   * @param point Label of point.
   *
   * @returns True if point is a constraint point, false otherwise.
   */
  bool _isConstraintEdge(const int point) const;

}; // class TestFaultCohesiveImpulses

#endif // pylith_faults_testfaultcohesiveimpulses_hh

// End of file 
