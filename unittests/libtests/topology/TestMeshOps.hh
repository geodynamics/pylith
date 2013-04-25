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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/topology/TestMeshOps.hh
 *
 * @brief C++ TestMeshOps object.
 * 
 * C++ unit testing for MeshOps.
 */

#if !defined(pylith_topology_testmeshops_hh)
#define pylith_topology_testmeshops_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestMeshOps;
  } // topology
} // pylith

/// C++ unit testing for MeshOps.
class pylith::topology::TestMeshOps : public CppUnit::TestFixture
{ // class TestMeshOps

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMeshOps );

  CPPUNIT_TEST( testCheckMaterialIds );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test checkMaterialIds().
  void testCheckMaterialIds(void);

}; // class TestMeshOps

#endif // pylith_topology_meshops_hh


// End of file 
