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
 * @file unittests/libtests/feassemble/TestElasticityExplicitLgDeformCases.hh
 *
 * @brief C++ TestElasticityExplicitLgDeform object
 *
 * C++ unit testing for ElasticityExplicitLgDeform with 1-D linear cells.
 */

#if !defined(pylith_feassemble_testelasticityexplicitlgdeform1dlinear_hh)
#define pylith_feassemble_testelasticityexplicitlgdeform1dlinear_hh

#include "TestElasticityExplicitLgDeform.hh" // ISA TestElasticityExplicitLgDeform

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityExplicitLgDeform2DLinear;
    class TestElasticityExplicitLgDeform2DQuadratic;
    class TestElasticityExplicitLgDeform3DLinear;
    class TestElasticityExplicitLgDeform3DQuadratic;

    class TestElasticityExplicitLgDeformGrav2DLinear;
    class TestElasticityExplicitLgDeformGrav2DQuadratic;
    class TestElasticityExplicitLgDeformGrav3DLinear;
    class TestElasticityExplicitLgDeformGrav3DQuadratic;
  } // feassemble
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicitLgDeform w/2-D linear cells.
class pylith::feassemble::TestElasticityExplicitLgDeform2DLinear :
  public TestElasticityExplicitLgDeform
{ // class TestElasticityExplicitLgDeform2DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitLgDeform2DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityExplicitLgDeform2DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicitLgDeform w/2-D quadratic cells.
class pylith::feassemble::TestElasticityExplicitLgDeform2DQuadratic :
  public TestElasticityExplicitLgDeform
{ // class TestElasticityExplicitLgDeform2DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitLgDeform2DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityExplicitLgDeform2DQuadratic


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicitLgDeform w/3-D linear cells.
class pylith::feassemble::TestElasticityExplicitLgDeform3DLinear :
  public TestElasticityExplicitLgDeform
{ // class TestElasticityExplicitLgDeform3DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitLgDeform3DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityExplicitLgDeform3DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicitLgDeform w/3-D quadratic cells.
class pylith::feassemble::TestElasticityExplicitLgDeform3DQuadratic :
  public TestElasticityExplicitLgDeform
{ // class TestElasticityExplicitLgDeform3DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitLgDeform3DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityExplicitLgDeform3DQuadratic


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicitLgDeform w/2-D linear cells and gravity.
class pylith::feassemble::TestElasticityExplicitLgDeformGrav2DLinear :
  public TestElasticityExplicitLgDeform
{ // class TestElasticityExplicitLgDeformGrav2DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitLgDeformGrav2DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityExplicitLgDeformGrav2DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicitLgDeform w/2-D quadratic cells and gravity.
class pylith::feassemble::TestElasticityExplicitLgDeformGrav2DQuadratic :
  public TestElasticityExplicitLgDeform
{ // class TestElasticityExplicitLgDeformGrav2DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitLgDeformGrav2DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityExplicitLgDeformGrav2DQuadratic


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicitLgDeform w/3-D linear cells and gravity.
class pylith::feassemble::TestElasticityExplicitLgDeformGrav3DLinear :
  public TestElasticityExplicitLgDeform
{ // class TestElasticityExplicitLgDeformGrav3DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitLgDeformGrav3DLinear );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityExplicitLgDeformGrav3DLinear


// ----------------------------------------------------------------------
/// C++ unit testing for ElasticityExplicitLgDeform w/3-D quadratic cells and gravity.
class pylith::feassemble::TestElasticityExplicitLgDeformGrav3DQuadratic :
  public TestElasticityExplicitLgDeform
{ // class TestElasticityExplicitLgDeformGrav3DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitLgDeformGrav3DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityExplicitLgDeformGrav3DQuadratic


#endif // pylith_feassemble_testelasticityexplicitlgdeformcases_hh


// End of file 
