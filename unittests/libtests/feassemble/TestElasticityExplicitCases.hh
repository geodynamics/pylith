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
 * @file unittests/libtests/feassemble/TestElasticityExplicitCases.hh
 *
 * @brief C++ TestElasticityExplicit object
 *
 * C++ unit testing for ElasticityExplicit.
 */

#if !defined(pylith_feassemble_testelasticityexplicitcases_hh)
#define pylith_feassemble_testelasticityexplicitcases_hh

#include "TestElasticityExplicit.hh" // ISA TestElasticityExplicit

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityExplicit2DLinear;
    class TestElasticityExplicit2DQuadratic;
    class TestElasticityExplicit3DLinear;
    class TestElasticityExplicit3DQuadratic;

    class TestElasticityExplicitGrav2DLinear;
    class TestElasticityExplicitGrav2DQuadratic;
    class TestElasticityExplicitGrav3DLinear;
    class TestElasticityExplicitGrav3DQuadratic;
  } // feassemble
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicit w/2-D linear cells.
class pylith::feassemble::TestElasticityExplicit2DLinear :
  public TestElasticityExplicit
{ // class TestElasticityExplicit2DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicit2DLinear );

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

}; // class TestElasticityExplicit2DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicit w/2-D quadratic cells.
class pylith::feassemble::TestElasticityExplicit2DQuadratic :
  public TestElasticityExplicit
{ // class TestElasticityExplicit2DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicit2DQuadratic );

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

}; // class TestElasticityExplicit2DQuadratic


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicit w/3-D linear cells.
class pylith::feassemble::TestElasticityExplicit3DLinear :
  public TestElasticityExplicit
{ // class TestElasticityExplicit3DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicit3DLinear );

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

}; // class TestElasticityExplicit3DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicit w/3-D quadratic cells.
class pylith::feassemble::TestElasticityExplicit3DQuadratic :
  public TestElasticityExplicit
{ // class TestElasticityExplicit3DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicit3DQuadratic );

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

}; // class TestElasticityExplicit3DQuadratic



// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicit w/2-D linear cells and gravity.
class pylith::feassemble::TestElasticityExplicitGrav2DLinear :
  public TestElasticityExplicit
{ // class TestElasticityExplicitGrav2DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitGrav2DLinear );

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

}; // class TestElasticityExplicitGrav2DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicit w/2-D quadratic cells amd gravity.
class pylith::feassemble::TestElasticityExplicitGrav2DQuadratic :
  public TestElasticityExplicit
{ // class TestElasticityExplicitGrav2DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitGrav2DQuadratic );

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

}; // class TestElasticityExplicitGrav2DQuadratic


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicit w/3-D linear cells and gravity.
class pylith::feassemble::TestElasticityExplicitGrav3DLinear :
  public TestElasticityExplicit
{ // class TestElasticityExplicitGrav3DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitGrav3DLinear );

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

}; // class TestElasticityExplicitGrav3DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicit w/3-D quadratic cells amd gravity.
class pylith::feassemble::TestElasticityExplicitGrav3DQuadratic :
  public TestElasticityExplicit
{ // class TestElasticityExplicitGrav3DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitGrav3DQuadratic );

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

}; // class TestElasticityExplicitGrav3DQuadratic


#endif // pylith_feassemble_testelasticityexplicitcases_hh


// End of file 
