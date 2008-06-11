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
 * @file unittests/libtests/feassemble/TestElasticityImplicitGrav2DQuadratic.hh
 *
 * @brief C++ TestElasticityImplicit object
 *
 * C++ unit testing for ElasticityImplicit with 2-D quadratic cells and gravity.
 */

#if !defined(pylith_feassemble_testelasticityimplicitgrav2dquadratic_hh)
#define pylith_feassemble_testelasticityimplicitgrav2dquadratic_hh

#include "TestElasticityImplicit.hh" // ISA TestElasticityImplicit

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityImplicitGrav2DQuadratic;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityImplicit
class pylith::feassemble::TestElasticityImplicitGrav2DQuadratic :
  public TestElasticityImplicit
{ // class TestElasticityImplicit2DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitGrav2DQuadratic );

  CPPUNIT_TEST( testUpdateState );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestElasticityImplicitGrav2DQuadratic

#endif // pylith_feassemble_testelasticityimplicitgrav2dquadratic_hh


// End of file 
