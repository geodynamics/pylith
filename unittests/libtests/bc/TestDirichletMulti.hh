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
 * @file unittests/libtests/bc/TestDirichletMulti.hh
 *
 * @brief C++ TestDirichletMulti object.
 *
 * C++ unit testing for DirichletMulti.
 */

#if !defined(pylith_bc_testdirichletmulti_hh)
#define pylith_bc_testdirichletmulti_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletMulti;

    class Dirichlet;
    class DirichletDataMulti;
  } // bc
} // pylith

/// C++ unit testing for DirichletMulti.
class pylith::bc::TestDirichletMulti : public CppUnit::TestFixture
{ // class TestDirichletMulti

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

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

  /** Initialize DirichletMulti boundary condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param bcA Dirichlet boundary condition A to initialize.
   * @param bcB Dirichlet boundary condition B to initialize.
   */
  void _initialize(ALE::Obj<ALE::Mesh>* mesh,
		   Dirichlet* const bcA,
		   Dirichlet* const bcB) const;

}; // class TestDirichletMulti

#endif // pylith_bc_dirichletmulti_hh


// End of file 
