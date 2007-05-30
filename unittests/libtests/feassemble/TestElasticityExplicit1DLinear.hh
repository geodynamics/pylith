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
 * @file unittests/libtests/feassemble/TestElasticityExplicit1DLinear.hh
 *
 * @brief C++ TestElasticityExplicit object
 *
 * C++ unit testing for ElasticityExplicit with 1-D linear cells.
 */

#if !defined(pylith_feassemble_testelasticityexplicit1dlinear_hh)
#define pylith_feassemble_testelasticityexplicit1dlinear_hh

#include "TestElasticityExplicit.hh" // ISA TestElasticityExplicit

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityExplicit1DLinear;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityExplicit
class pylith::feassemble::TestElasticityExplicit1DLinear :
  public TestElasticityExplicit
{ // class TestElasticityExplicit1DLinear

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicit1DLinear );

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

}; // class TestElasticityExplicit1DLinear

#endif // pylith_feassemble_testelasticityexplicit1dlinear_hh


// End of file 
