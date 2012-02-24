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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/feassemble/TestElasticityExplicit1DLinear.hh
 *
 * @brief C++ TestElasticityExplicit object
 *
 * C++ unit testing for ElasticityExplicit with 1-D linear cells.
 */

#if !defined(pylith_feassemble_testelasticityexplicit1dlinear_hh)
#define pylith_feassemble_testelasticityexplicit1dlinear_hh

#include "TestElasticityExplicit.hh" // ISA TestElasticityExplicit

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityExplicit1DLinear;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityExplicit
class pylith::feassemble::TestElasticityExplicit1DLinear :
  public TestElasticityExplicit
{ // class TestElasticityExplicit1DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicit1DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateResidualLumped );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testStableTimeStep );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityExplicit1DLinear

#endif // pylith_feassemble_testelasticityexplicit1dlinear_hh


// End of file 
