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
 * @file unittests/libtests/feassemble/TestIntegratorElasticity.hh
 *
 * @brief C++ TestIntegratorElasticity object
 *
 * C++ unit testing for IntegratorElasticity.
 */

#if !defined(pylith_feassemble_testintegratorelasticity_hh)
#define pylith_feassemble_testintegratorelasticity_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestIntegratorElasticity;
  } // feassemble
} // pylith

/// C++ unit testing for Elasticity
class pylith::feassemble::TestIntegratorElasticity : public CppUnit::TestFixture
{ // class TestIntegratorElasticity

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIntegratorElasticity );

  CPPUNIT_TEST( testCalcTotalStrain2D );
  CPPUNIT_TEST( testCalcTotalStrain3D );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test calcTotalStrain2D().
  void testCalcTotalStrain2D(void);

  /// Test calcTotalStrain3D().
  void testCalcTotalStrain3D(void);

}; // class TestIntegratorElasticity

#endif // pylith_feassemble_testintegratorelasticity_hh


// End of file 
