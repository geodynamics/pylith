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
 * @file unittests/libtests/topology/TestSolutionFields.hh
 *
 * @brief C++ TestSolutionFields object.
 * 
 * C++ unit testing for SolutionFields.
 */

#if !defined(pylith_topology_testsolutionfields_hh)
#define pylith_topology_testsolutionfields_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh"

/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestSolutionFields;
  } // topology
} // pylith

/// C++ unit testing for SolutionFields.
class pylith::topology::TestSolutionFields : public CppUnit::TestFixture
{ // class TestSolutionFields

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestSolutionFields );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testSolutionName );
  CPPUNIT_TEST( testSolution );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test solutionName().
  void testSolutionName(void);

  /// Test solution().
  void testSolution(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize mesh for SolutionFields.
   *
   * @param mesh Finite-element mesh.
   */
  void _initialize(Mesh* mesh) const;

}; // class TestSolutionFields

#endif // pylith_topology_solutionfields_hh


// End of file 
