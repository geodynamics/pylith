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
 * @file unittests/libtests/bc/TestAbsorbingDampers.hh
 *
 * @brief C++ TestAbsorbingDampers object.
 *
 * C++ unit testing for AbsorbingDampers.
 */

#if !defined(pylith_bc_testabsorbingdampers_hh)
#define pylith_bc_testabsorbingdampers_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/bc/bcfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/feassemble/feassemblefwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestAbsorbingDampers;
    class AbsorbingDampersData;
  } // bc
} // pylith

/// C++ unit testing for AbsorbingDampers.
class pylith::bc::TestAbsorbingDampers : public CppUnit::TestFixture
{ // class TestAbsorbingDampers

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestAbsorbingDampers );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testDB );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Test constructor.
  void testConstructor(void);

  /// Test db()
  void testDB(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test integrateResidual().
  void testIntegrateResidual(void);

  /// Test integrateJacobian().
  void testIntegrateJacobian(void);

  /// Test integrateJacobianLumped().
  void testIntegrateJacobianLumped(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  AbsorbingDampersData* _data; ///< Data for testing
  feassemble::Quadrature* _quadrature; ///< Used in testing.

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize AbsorbingDampers boundary condition.
   *
   * @param mesh Finite-element mesh to initialize
   * @param bc Neumann boundary condition to initialize.
   * @param fields Solution fields.
   */
  void _initialize(topology::Mesh* mesh,
		   AbsorbingDampers* const bc,
		   topology::SolutionFields* fields) const;

}; // class TestAbsorbingDampers

#endif // pylith_bc_absorbingdampers_hh


// End of file 
