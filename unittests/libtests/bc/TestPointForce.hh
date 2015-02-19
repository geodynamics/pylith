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
 * @file unittests/libtests/bc/TestPointForce.hh
 *
 * @brief C++ TestPointForce object.
 *
 * C++ unit testing for PointForce.
 */

#if !defined(pylith_bc_testpointforce_hh)
#define pylith_bc_testpointforce_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/bc/bcfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestPointForce;
    class PointForceData;
  } // bc
} // pylith

/// C++ unit testing for PointForce.
class pylith::bc::TestPointForce : public CppUnit::TestFixture
{ // class TestPointForce

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestPointForce );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testNormalizer );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Test constructor.
  void testConstructor(void);

  /// Test normalizer().
  void testNormalizer(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test integrateResidual().
  void testIntegrateResidual(void);

  /// Test verifyConfiguration().
  void testVerifyConfiguration(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  PointForceData* _data; ///< Data for testing

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize PointForce boundary condition.
   *
   * @param mesh Finite-element mesh to initialize.
   * @param bc PointForce boundary condition to initialize.
   */
  void _initialize(topology::Mesh* mesh,
		   PointForce* const bc) const;

}; // class TestPointForce

#endif // pylith_bc_pointforce_hh


// End of file 
