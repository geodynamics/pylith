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
 * @file unittests/libtests/bc/TestDirichletBoundary.hh
 *
 * @brief C++ TestDirichletBoundary object.
 *
 * C++ unit testing for DirichletBoundary.
 */

#if !defined(pylith_bc_testdirichletboundary_hh)
#define pylith_bc_testdirichletboundary_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/bc/bcfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBoundary;
    class DirichletData;
  } // bc

  namespace topology {
    class Mesh; // USES Mesh
  } // topology
} // pylith

/// C++ unit testing for DirichletBoundary.
class pylith::bc::TestDirichletBoundary : public CppUnit::TestFixture
{ // class TestDirichletBoundary

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletBoundary );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Test constructor.
  void testConstructor(void);

  /// Test initialize().
  void testInitialize(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  DirichletData* _data; ///< Data for testing

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize DirichletBoundary boundary condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param bc DirichletBoundary boundary condition to initialize.
   */
  void _initialize(topology::Mesh* mesh,
		   DirichletBoundary* const bc) const;

}; // class TestDirichletBoundary

#endif // pylith_bc_dirichletboundary_hh


// End of file 
