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
 * @file unittests/libtests/feassemble/TestElasticityExplicitLgDeformGrav2DLinear.hh
 *
 * @brief C++ TestElasticityExplicitLgDeform object
 *
 * C++ unit testing for ElasticityExplicitLgDeform with 2-D linear cells and gravity.
 */

#if !defined(pylith_feassemble_testelasticityexplicitlgdeformgrav2dlinear_hh)
#define pylith_feassemble_testelasticityexplicitlgdeformgrav2dlinear_hh

#include "TestElasticityExplicitLgDeform.hh" // ISA TestElasticityExplicitLgDeform

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityExplicitLgDeformGrav2DLinear;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityExplicitLgDeform
class pylith::feassemble::TestElasticityExplicitLgDeformGrav2DLinear :
  public TestElasticityExplicitLgDeform
{ // class TestElasticityExplicitLgDeform2DLinear

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

#endif // pylith_feassemble_testelasticityexplicitlgdeformgrav2dlinear_hh


// End of file 
