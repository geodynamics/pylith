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
 * @file unittests/libtests/bc/TestDirichletBoundary.hh
 *
 * @brief C++ TestDirichletBoundary object.
 *
 * C++ unit testing for DirichletBoundary.
 */

#if !defined(pylith_bc_testdirichletboundary_hh)
#define pylith_bc_testdirichletboundary_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBoundary;

    class DirichletBoundary;
    class DirichletData;
  } // bc
} // pylith

/// C++ unit testing for DirichletBoundary.
class pylith::bc::TestDirichletBoundary : public CppUnit::TestFixture
{ // class TestDirichletBoundary

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletBoundary );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testFixedDOF );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

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

  /** Initialize DirichletBoundary boundary condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param bc DirichletBoundary boundary condition to initialize.
   */
  void _initialize(ALE::Obj<ALE::Mesh>* mesh,
		   DirichletBoundary* const bc) const;

}; // class TestDirichletBoundary

#endif // pylith_bc_dirichletboundary_hh


// End of file 
