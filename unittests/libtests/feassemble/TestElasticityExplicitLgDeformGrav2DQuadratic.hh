// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/feassemble/TestElasticityExplicitLgDeformGrav2DQuadratic.hh
 *
 * @brief C++ TestElasticityExplicitLgDeform object
 *
 * C++ unit testing for ElasticityExplicitLgDeform with 2-D quadratic cells and gravity.
 */

#if !defined(pylith_feassemble_testelasticityexplicitlgdeformgrav2dquadratic_hh)
#define pylith_feassemble_testelasticityexplicitlgdeformgrav2dquadratic_hh

#include "TestElasticityExplicitLgDeform.hh" // ISA TestElasticityExplicitLgDeform

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityExplicitLgDeformGrav2DQuadratic;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityExplicitLgDeform
class pylith::feassemble::TestElasticityExplicitLgDeformGrav2DQuadratic :
  public TestElasticityExplicitLgDeform
{ // class TestElasticityExplicitLgDeform2DQuadratic

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

#endif // pylith_feassemble_testelasticityexplicitlgdeformgrav2dquadratic_hh


// End of file 
