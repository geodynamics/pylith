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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/feassemble/TestElasticityImplicitGrav3DQuadratic.hh
 *
 * @brief C++ TestElasticityImplicit object
 *
 * C++ unit testing for ElasticityImplicit with 1-D quadratic cells and gravity.
 */

#if !defined(pylith_feassemble_testelasticityimplicitgrav3dquadratic_hh)
#define pylith_feassemble_testelasticityimplicitgrav3dquadratic_hh

#include "TestElasticityImplicit.hh" // ISA TestElasticityImplicit

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityImplicitGrav3DQuadratic;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityImplicit
class pylith::feassemble::TestElasticityImplicitGrav3DQuadratic :
  public TestElasticityImplicit
{ // class TestElasticityImplicit3DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitGrav3DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitGrav3DQuadratic

#endif // pylith_feassemble_testelasticityimplicitgrav3dquadratic_hh


// End of file 
