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
 * @file unittests/libtests/bc/TestDirichletBC.hh
 *
 * @brief C++ TestDirichletBC object.
 *
 * C++ unit testing for DirichletBC.
 */

#if !defined(pylith_bc_testdirichletbc_hh)
#define pylith_bc_testdirichletbc_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBC;

    class DirichletBC;
    class DirichletData;
  } // bc

  namespace topology {
    class Mesh;
  } // topology
} // pylith

/// C++ unit testing for DirichletBC.
class pylith::bc::TestDirichletBC : public CppUnit::TestFixture
{ // class TestDirichletBC

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletBC );
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
