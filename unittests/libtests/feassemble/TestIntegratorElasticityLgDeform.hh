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
 * @file unittests/libtests/feassemble/TestIntegratorElasticityLgDeform.hh
 *
 * @brief C++ TestIntegratorElasticityLgDeform object
 *
 * C++ unit testing for IntegratorElasticityLgDeform.
 */

#if !defined(pylith_feassemble_testintegratorelasticitylgdeform_hh)
#define pylith_feassemble_testintegratorelasticitylgdeform_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestIntegratorElasticityLgDeform;
  } // feassemble
} // pylith

/// C++ unit testing for Elasticity
class pylith::feassemble::TestIntegratorElasticityLgDeform : public CppUnit::TestFixture
{ // class TestIntegratorElasticityLgDeform

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIntegratorElasticityLgDeform );

  CPPUNIT_TEST( testCalcDeformation2D );
  CPPUNIT_TEST( testCalcDeformation3D );
  CPPUNIT_TEST( testCalcTotalStrain2D );
  CPPUNIT_TEST( testCalcTotalStrain3D );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test calcDeformation() for 2-D.
  void testCalcDeformation2D(void);

  /// Test calcDeformation() for 3-D.
  void testCalcDeformation3D(void);

  /// Test calcTotalStrain2D().
  void testCalcTotalStrain2D(void);

  /// Test calcTotalStrain3D().
  void testCalcTotalStrain3D(void);

}; // class TestIntegratorElasticityLgDeform

#endif // pylith_feassemble_testintegratorelasticitylgdeform_hh


// End of file 
