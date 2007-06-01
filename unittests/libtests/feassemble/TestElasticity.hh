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
 * @file unittests/libtests/feassemble/TestElasticity.hh
 *
 * @brief C++ TestElasticity object
 *
 * C++ unit testing for Elasticity.
 */

#if !defined(pylith_feassemble_testelasticity_hh)
#define pylith_feassemble_testelasticity_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticity;
  } // feassemble
} // pylith

/// C++ unit testing for Elasticity
class pylith::feassemble::TestElasticity : public CppUnit::TestFixture
{ // class TestElasticity

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticity );

  CPPUNIT_TEST( testCalcTotalStrain1D );
  CPPUNIT_TEST( testCalcTotalStrain2D );
  CPPUNIT_TEST( testCalcTotalStrain3D );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test calcTotalStrain1D().
  void testCalcTotalStrain1D(void);

  /// Test calcTotalStrain2D().
  void testCalcTotalStrain2D(void);

  /// Test calcTotalStrain3D().
  void testCalcTotalStrain3D(void);

}; // class TestElasticity

#endif // pylith_feassemble_testelasticity_hh


// End of file 
