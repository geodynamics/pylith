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
 * @file unittests/libtests/bc/TestDirichletBoundaryMulti.hh
 *
 * @brief C++ TestDirichletBoundaryMulti object.
 *
 * C++ unit testing for DirichletBoundaryMulti.
 */

#if !defined(pylith_bc_testdirichletboundarymulti_hh)
#define pylith_bc_testdirichletboundarymulti_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBoundaryMulti;

    class DirichletBoundary;
    class DirichletDataMulti;
  } // bc
} // pylith

/// C++ unit testing for DirichletBoundaryMulti.
class pylith::bc::TestDirichletBoundaryMulti : public CppUnit::TestFixture
{ // class TestDirichletBoundaryMulti

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Test setConstraintSizes().
  void testSetConstraintSizes(void);

  /// Test setConstraints().
  void testSetConstraints(void);

  /// Test setField().
  void testSetField(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  DirichletDataMulti* _data; ///< Data for testing

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize DirichletBoundaryMulti boundary condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param bcA DirichletBoundary boundary condition A to initialize.
   * @param bcB DirichletBoundary boundary condition B to initialize.
   */
  void _initialize(ALE::Obj<ALE::Mesh>* mesh,
		   DirichletBoundary* const bcA,
		   DirichletBoundary* const bcB) const;

}; // class TestDirichletBoundaryMulti

#endif // pylith_bc_dirichletboundarymulti_hh


// End of file 
