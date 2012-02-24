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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/bc/TestAbsorbingDampersTri3.hh
 *
 * @brief C++ TestAbsorbingDampers object.
 *
 * C++ unit testing for AbsorbingDampers for mesh with 2-D tri cells.
 */

#if !defined(pylith_bc_testabsorbingdamperstri3_hh)
#define pylith_bc_testabsorbingdamperstri3_hh

#include "TestAbsorbingDampers.hh" // ISA TestAbsorbingDampers

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestAbsorbingDampersTri3;
  } // bc
} // pylith

/// C++ unit testing for AbsorbingDampers for mesh with 2-D tri cells.
class pylith::bc::TestAbsorbingDampersTri3 : public TestAbsorbingDampers
{ // class TestAbsorbingDampers

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestAbsorbingDampersTri3, TestAbsorbingDampers );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestAbsorbingDampersTri3

#endif // pylith_bc_absorbingdamperstri3_hh


// End of file 
