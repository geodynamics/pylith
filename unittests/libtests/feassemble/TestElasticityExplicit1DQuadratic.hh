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
 * @file unittests/libtests/feassemble/TestElasticityExplicit1DQuadratic.hh
 *
 * @brief C++ TestElasticityExplicit object
 *
 * C++ unit testing for ElasticityExplicit with 1-D quadratic cells.
 */

#if !defined(pylith_feassemble_testelasticityexplicit1dquadratic_hh)
#define pylith_feassemble_testelasticityexplicit1dquadratic_hh

#include "TestElasticityExplicit.hh" // ISA TestElasticityExplicit

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityExplicit1DQuadratic;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityExplicit
class pylith::feassemble::TestElasticityExplicit1DQuadratic :
  public TestElasticityExplicit
{ // class TestElasticityExplicit1DQuadratic

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicit1DQuadratic );

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

}; // class TestElasticityExplicit1DQuadratic

#endif // pylith_feassemble_testelasticityexplicit1dquadratic_hh


// End of file 
