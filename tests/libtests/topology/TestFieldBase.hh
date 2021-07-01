// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file tests/libtests/topology/TestFieldBase.hh
 *
 * @brief C++ unit testing for Field.
 */

#if !defined(pylith_topology_testfieldbase_hh)
#define pylith_topology_testfieldbase_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // forward declarations

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestFieldBase;
  } // topology
} // pylith

// TestFieldBase -------------------------------------------------------------
/// C++ unit testing for Field.
class pylith::topology::TestFieldBase : public CppUnit::TestFixture
{ // class TestFieldBase

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFieldBase );

  CPPUNIT_TEST( testVectorFieldString );
  CPPUNIT_TEST( testParseVectorFieldString );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test vectorFieldString().
  void testVectorFieldString(void);

  /// Test parseVectorFieldString().
  void testParseVectorFieldString(void);


}; // class TestFieldBase

#endif // pylith_topology_testfieldbase_hh


// End of file 
