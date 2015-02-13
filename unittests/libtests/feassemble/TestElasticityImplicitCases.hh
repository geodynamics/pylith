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
 * @file unittests/libtests/feassemble/TestElasticityImplicitCases.hh
 *
 * @brief C++ TestElasticityImplicit object
 *
 * C++ unit testing for ElasticityImplicit.
 */

#if !defined(pylith_feassemble_testelasticityimplicitcases_hh)
#define pylith_feassemble_testelasticityimplicitcases_hh

#include "TestElasticityImplicit.hh" // ISA TestElasticityImplicit

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityImplicit2DLinear;
    class TestElasticityImplicit2DQuadratic;
    class TestElasticityImplicit3DLinear;
    class TestElasticityImplicit3DQuadratic;

    class TestElasticityImplicitGrav2DLinear;
    class TestElasticityImplicitGrav2DQuadratic;
    class TestElasticityImplicitGrav3DLinear;
    class TestElasticityImplicitGrav3DQuadratic;
  } // feassemble
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicit w/2-D linear cells.
class pylith::feassemble::TestElasticityImplicit2DLinear :
  public TestElasticityImplicit
{ // class TestElasticityImplicit2DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicit2DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testStableTimeStep );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicit2DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicit w/2-D quadratic cells.
class pylith::feassemble::TestElasticityImplicit2DQuadratic :
  public TestElasticityImplicit
{ // class TestElasticityImplicit2DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicit2DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testStableTimeStep );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicit2DQuadratic


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicit w/3-D linear cells.
class pylith::feassemble::TestElasticityImplicit3DLinear :
  public TestElasticityImplicit
{ // class TestElasticityImplicit3DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicit3DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testStableTimeStep );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicit3DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicit w/3-D quadratic cells.
class pylith::feassemble::TestElasticityImplicit3DQuadratic :
  public TestElasticityImplicit
{ // class TestElasticityImplicit3DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicit3DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testStableTimeStep );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicit3DQuadratic


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicit w/2-D linear cells and gravity.
class pylith::feassemble::TestElasticityImplicitGrav2DLinear :
  public TestElasticityImplicit
{ // class TestElasticityImplicitGrav2DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitGrav2DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testStableTimeStep );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitGrav2DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicit w/2-D quadratic cells and gravity.
class pylith::feassemble::TestElasticityImplicitGrav2DQuadratic :
  public TestElasticityImplicit
{ // class TestElasticityImplicitGrav2DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitGrav2DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testStableTimeStep );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitGrav2DQuadratic


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicit w/3-D linear cells and gravity.
class pylith::feassemble::TestElasticityImplicitGrav3DLinear :
  public TestElasticityImplicit
{ // class TestElasticityImplicitGrav3DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitGrav3DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testStableTimeStep );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitGrav3DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicit w/3-D quadratic cells and gravity.
class pylith::feassemble::TestElasticityImplicitGrav3DQuadratic :
  public TestElasticityImplicit
{ // class TestElasticityImplicitGrav3DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitGrav3DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testStableTimeStep );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitGrav3DQuadratic


#endif // pylith_feassemble_testelasticityimplicitcases_hh


// End of file 
