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
 * @file unittests/libtests/feassemble/TestElasticityExplicitLgDeformGrav1DQuadratic.hh
 *
 * @brief C++ TestElasticityExplicitLgDeform object
 *
 * C++ unit testing for ElasticityExplicitLgDeform with 1-D quadratic cells and gravity.
 */

#if !defined(pylith_feassemble_testelasticityexplicitlgdeformgrav1dquadratic_hh)
#define pylith_feassemble_testelasticityexplicitlgdeformgrav1dquadratic_hh

#include "TestElasticityExplicitLgDeform.hh" // ISA TestElasticityExplicitLgDeform

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityExplicitLgDeformGrav1DQuadratic;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityExplicitLgDeform
class pylith::feassemble::TestElasticityExplicitLgDeformGrav1DQuadratic :
  public TestElasticityExplicitLgDeform
{ // class TestElasticityExplicitLgDeform1DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitLgDeformGrav1DQuadratic );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateResidualLumped );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityExplicitLgDeformGrav1DQuadratic

#endif // pylith_feassemble_testelasticityexplicitlgdeformgrav1dquadratic_hh


// End of file 
