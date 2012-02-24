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
 * @file unittests/libtests/feassemble/TestElasticityImplicitLgDeform1DQuadratic.hh
 *
 * @brief C++ TestElasticityImplicitLgDeform object
 *
 * C++ unit testing for ElasticityImplicitLgDeform with 1-D quadratic cells.
 */

#if !defined(pylith_feassemble_testelasticityimplicitlgdeform1dquadratic_hh)
#define pylith_feassemble_testelasticityimplicitlgdeform1dquadratic_hh

#include "TestElasticityImplicitLgDeform.hh" // ISA TestElasticityImplicitLgDeform

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityImplicitLgDeform1DQuadratic;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityImplicitLgDeform
class pylith::feassemble::TestElasticityImplicitLgDeform1DQuadratic :
  public TestElasticityImplicitLgDeform
{ // class TestElasticityImplicitLgDeform1DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitLgDeform1DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitLgDeform1DQuadratic

#endif // pylith_feassemble_testelasticityimplicitlgdeform1dquadratic_hh


// End of file 
