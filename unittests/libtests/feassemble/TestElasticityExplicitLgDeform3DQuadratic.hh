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
 * @file unittests/libtests/feassemble/TestElasticityExplicitLgDeform3DQuadratic.hh
 *
 * @brief C++ TestElasticityExplicitLgDeform object
 *
 * C++ unit testing for ElasticityExplicitLgDeform with 1-D quadratic cells.
 */

#if !defined(pylith_feassemble_testelasticityexplicitlgdeform3dquadratic_hh)
#define pylith_feassemble_testelasticityexplicitlgdeform3dquadratic_hh

#include "TestElasticityExplicitLgDeform.hh" // ISA TestElasticityExplicitLgDeform

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityExplicitLgDeform3DQuadratic;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityExplicitLgDeform
class pylith::feassemble::TestElasticityExplicitLgDeform3DQuadratic :
  public TestElasticityExplicitLgDeform
{ // class TestElasticityExplicitLgDeform3DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitLgDeform3DQuadratic );

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

}; // class TestElasticityExplicitLgDeform3DQuadratic

#endif // pylith_feassemble_testelasticityexplicitlgdeform3dquadratic_hh


// End of file 
