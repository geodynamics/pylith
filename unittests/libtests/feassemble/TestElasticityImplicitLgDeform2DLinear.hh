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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/feassemble/TestElasticityImplicitLgDeform2DLinear.hh
 *
 * @brief C++ TestElasticityImplicitLgDeform object
 *
 * C++ unit testing for ElasticityImplicitLgDeform with 2-D linear cells.
 */

#if !defined(pylith_feassemble_testelasticityimplicitlgdeform2dlinear_hh)
#define pylith_feassemble_testelasticityimplicitlgdeform2dlinear_hh

#include "TestElasticityImplicitLgDeform.hh" // ISA TestElasticityImplicitLgDeform

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityImplicitLgDeform2DLinear;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityImplicitLgDeform
class pylith::feassemble::TestElasticityImplicitLgDeform2DLinear :
  public TestElasticityImplicitLgDeform
{ // class TestElasticityImplicitLgDeform2DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitLgDeform2DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitLgDeform2DLinear

#endif // pylith_feassemble_testelasticityimplicitlgdeform2dlinear_hh


// End of file 
