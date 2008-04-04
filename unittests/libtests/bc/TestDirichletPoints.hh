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
 * @file unittests/libtests/bc/TestDirichletPoints.hh
 *
 * @brief C++ TestDirichletPoints object.
 *
 * C++ unit testing for DirichletPoints.
 */

#if !defined(pylith_bc_testdirichletpoints_hh)
#define pylith_bc_testdirichletpoints_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletPoints;

    class DirichletPoints;
    class DirichletData;
  } // bc
} // pylith

/// C++ unit testing for DirichletPoints.
class pylith::bc::TestDirichletPoints : public CppUnit::TestFixture
{ // class TestDirichletPoints

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletPoints );
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

  /** Initialize DirichletPoints boundary condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param bc DirichletPoints boundary condition to initialize.
   */
  void _initialize(ALE::Obj<Mesh>* mesh,
		   DirichletPoints* const bc) const;

}; // class TestDirichletPoints

#endif // pylith_bc_dirichletpoints_hh


// End of file 
