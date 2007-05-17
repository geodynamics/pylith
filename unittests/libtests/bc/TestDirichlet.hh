// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/bc/TestDirichlet.hh
 *
 * @brief C++ TestDirichlet object.
 *
 * C++ unit testing for Dirichlet.
 */

#if !defined(pylith_bc_testdirichlet_hh)
#define pylith_bc_testdirichlet_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichlet;

    class Dirichlet;
    class DirichletData;
  } // bc
} // pylith

/// C++ unit testing for Dirichlet.
class pylith::bc::TestDirichlet : public CppUnit::TestFixture
{ // class TestDirichlet

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichlet );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testFixedDOF );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test fixedDOF().
  void testFixedDOF(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test setConstraintSizes().
  void testSetConstraintSizes(void);

  /// Test setConstraints().
  void testSetConstraints(void);

  /// Test setField().
  void testSetField(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  DirichletData* _data; ///< Data for testing

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize Dirichlet boundary condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param bc Dirichlet boundary condition to initialize.
   */
  void _initialize(ALE::Obj<ALE::Mesh>* mesh,
		   Dirichlet* const bc) const;

}; // class TestDirichlet

#endif // pylith_bc_dirichlet_hh


// End of file 
