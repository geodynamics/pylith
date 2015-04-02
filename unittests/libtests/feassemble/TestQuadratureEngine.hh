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
 * @file unittests/libtests/feassemble/TestQuadratureEngine.hh
 *
 * @brief C++ TestQuadratureEngine object
 *
 * C++ unit testing for Quadrature.
 */

#if !defined(pylith_feassemble_testquadrature_hh)
#define pylith_feassemble_testquadrature_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/feassemble/feassemblefwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestQuadratureEngine;
    class QuadratureData;
  } // feassemble
} // pylith

/// C++ unit testing for Quadrature
class pylith::feassemble::TestQuadratureEngine : public CppUnit::TestFixture
{ // class TestQuadratureEngine

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestQuadratureEngine );

  CPPUNIT_TEST( testCopyConstructor );
  CPPUNIT_TEST( testInitialize );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test copy constructor.
  void testCopyConstructor(void);

  /// Test initialize().
  void testInitialize(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Test computeGeometry() and retrieveGeometry().
   *
   * @param engine Quadrature engine.
   * @param refCell Quadrature reference cell information.
   * @param data Data for testing quadrature
   */
  void _testComputeGeometry(QuadratureEngine* engine,
			    QuadratureRefCell* refCell,
			    const QuadratureData& data) const;

}; // class TestQuadratureEngine

#endif // pylith_feassemble_testquadratureengine_hh

// End of file 
