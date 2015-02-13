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
 * @file unittests/libtests/topology/TestJacobian.hh
 *
 * @brief C++ TestJacobian object.
 * 
 * C++ unit testing for Jacobian.
 */

#if !defined(pylith_topology_testjacobian_hh)
#define pylith_topology_testjacobian_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh"

/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestJacobian;
  } // topology
} // pylith

/// C++ unit testing for Jacobian.
class pylith::topology::TestJacobian : public CppUnit::TestFixture
{ // class TestJacobian

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestJacobian );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testConstructorSubDomain );
  CPPUNIT_TEST( testMatrix );
  CPPUNIT_TEST( testAssemble );
  CPPUNIT_TEST( testZero );
  CPPUNIT_TEST( testView );
  CPPUNIT_TEST( testWrite );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test constructor with subdomain.
  void testConstructorSubDomain(void);

  /// Test matrix().
  void testMatrix(void);

  /// Test assemble().
  void testAssemble(void);

  /// Test zero().
  void testZero(void);

  /// Test view().
  void testView(void);

  /// Test write().
  void testWrite(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize mesh for Jacobian.
   *
   * @param mesh Finite-element mesh.
   * @param field Solution field.
   */
  void _initializeMesh(Mesh* mesh) const;

  /** Initialize field for Jacobian.
   *
   * @param mesh Finite-element mesh.
   * @param field Solution field.
   */
  void _initializeField(Mesh* mesh,
                        Field* field) const;

}; // class TestJacobian

#endif // pylith_topology_jacobian_hh


// End of file 
