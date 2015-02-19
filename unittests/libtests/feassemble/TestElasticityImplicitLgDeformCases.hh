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
 * @file unittests/libtests/feassemble/TestElasticityImplicitLgDeformCases.hh
 *
 * @brief C++ TestElasticityImplicitLgDeform object
 *
 * C++ unit testing for ElasticityImplicitLgDeform.
 */

#if !defined(pylith_feassemble_testelasticityimplicitlgdeformcases_hh)
#define pylith_feassemble_testelasticityimplicitlgdeformcases_hh

#include "TestElasticityImplicitLgDeform.hh" // ISA TestElasticityImplicitLgDeform

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityImplicitLgDeform2DLinear;
    class TestElasticityImplicitLgDeform2DQuadratic;
    class TestElasticityImplicitLgDeform3DLinear;
    class TestElasticityImplicitLgDeform3DQuadratic;

    class TestElasticityImplicitLgDeformGrav2DLinear;
    class TestElasticityImplicitLgDeformGrav2DQuadratic;
    class TestElasticityImplicitLgDeformGrav3DLinear;
    class TestElasticityImplicitLgDeformGrav3DQuadratic;
  } // feassemble
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicitLgDeform w/2-D linear cells.
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


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicitLgDeform w/2-D quadratic cells.
class pylith::feassemble::TestElasticityImplicitLgDeform2DQuadratic :
  public TestElasticityImplicitLgDeform
{ // class TestElasticityImplicitLgDeform2DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitLgDeform2DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitLgDeform2DQuadratic


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicitLgDeform w/3-D linear cells.
class pylith::feassemble::TestElasticityImplicitLgDeform3DLinear :
  public TestElasticityImplicitLgDeform
{ // class TestElasticityImplicitLgDeform3DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitLgDeform3DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitLgDeform3DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicitLgDeform w/3-D quadratic cells.
class pylith::feassemble::TestElasticityImplicitLgDeform3DQuadratic :
  public TestElasticityImplicitLgDeform
{ // class TestElasticityImplicitLgDeform3DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitLgDeform3DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitLgDeform3DQuadratic


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicitLgDeform w/2-D linear cells and gravity.
class pylith::feassemble::TestElasticityImplicitLgDeformGrav2DLinear :
  public TestElasticityImplicitLgDeform
{ // class TestElasticityImplicitLgDeformGrav2DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitLgDeformGrav2DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitLgDeformGrav2DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicitLgDeform w/2-D quadratic cells and gravity.
class pylith::feassemble::TestElasticityImplicitLgDeformGrav2DQuadratic :
  public TestElasticityImplicitLgDeform
{ // class TestElasticityImplicitLgDeformGrav2DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitLgDeformGrav2DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitLgDeformGrav2DQuadratic


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicitLgDeform w/3-D linear cells and gravity.
class pylith::feassemble::TestElasticityImplicitLgDeformGrav3DLinear :
  public TestElasticityImplicitLgDeform
{ // class TestElasticityImplicitLgDeformGrav3DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitLgDeformGrav3DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitLgDeformGrav3DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityImplicitLgDeform w/3-D quadratic cells and gravity.
class pylith::feassemble::TestElasticityImplicitLgDeformGrav3DQuadratic :
  public TestElasticityImplicitLgDeform
{ // class TestElasticityImplicitLgDeformGrav3DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitLgDeformGrav3DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitLgDeformGrav3DQuadratic


#endif // pylith_feassemble_testelasticityimplicitlgdeformcases_hh


// End of file 
