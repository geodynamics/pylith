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
 * @file unittests/libtests/feassemble/TestElasticityImplicit3DQuadratic.hh
 *
 * @brief C++ TestElasticityImplicit object
 *
 * C++ unit testing for ElasticityImplicit with 3-D quadratic cells.
 */

#if !defined(pylith_feassemble_testelasticityimplicit3dquadratic_hh)
#define pylith_feassemble_testelasticityimplicit3dquadratic_hh

#include "TestElasticityImplicit.hh" // ISA TestElasticityImplicit

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityImplicit3DQuadratic;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityImplicit
class pylith::feassemble::TestElasticityImplicit3DQuadratic :
  public TestElasticityImplicit
{ // class TestElasticityImplicit3DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicit3DQuadratic );

  CPPUNIT_TEST( testUpdateState );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

}; // class TestElasticityImplicit3DQuadratic

#endif // pylith_feassemble_testelasticityimplicit3dquadratic_hh


// End of file 
