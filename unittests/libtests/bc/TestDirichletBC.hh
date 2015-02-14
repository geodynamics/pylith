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
 * @file unittests/libtests/bc/TestDirichletBC.hh
 *
 * @brief C++ TestDirichletBC object.
 *
 * C++ unit testing for DirichletBC.
 */

#if !defined(pylith_bc_testdirichletbc_hh)
#define pylith_bc_testdirichletbc_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/bc/bcfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBC;
    class DirichletData;
  } // bc
} // pylith

/// C++ unit testing for DirichletBC.
class pylith::bc::TestDirichletBC : public CppUnit::TestFixture
{ // class TestDirichletBC

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletBC );

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

  /// Test numDimConstrained().
  void testNumDimConstrained(void);

  /// Test setConstraintSizes().
  void testSetConstraintSizes(void);

  /// Test setConstraints().
  void testSetConstraints(void);

  /// Test setField().
  void testSetField(void);

  /// Test setFieldIncr().
  void testSetFieldIncr(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  DirichletData* _data; ///< Data for testing

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize DirichletBC boundary condition.
   *
   * @param mesh Finite-element mesh to initialize.
   * @param bc DirichletBC boundary condition to initialize.
   */
  void _initialize(topology::Mesh* mesh,
		   DirichletBC* const bc) const;

}; // class TestDirichletBC

#endif // pylith_bc_dirichletbc_hh


// End of file 
