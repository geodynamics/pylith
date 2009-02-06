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
 * @file unittests/libtests/bc/TestDirichletBCMulti.hh
 *
 * @brief C++ TestDirichletBCMulti object.
 *
 * C++ unit testing for DirichletBCMulti.
 */

#if !defined(pylith_bc_testdirichletbcmulti_hh)
#define pylith_bc_testdirichletbcmulti_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/bc/bcfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBCMulti;
    class DirichletDataMulti;
  } // bc
} // pylith

/// C++ unit testing for DirichletBCMulti.
class pylith::bc::TestDirichletBCMulti : public CppUnit::TestFixture
{ // class TestDirichletBCMulti

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

  /** Initialize DirichletBCMulti boundary condition.
   *
   * @param mesh Finite-element mesh to initialize.
   * @param bcA DirichletBC boundary condition A to initialize.
   * @param bcB DirichletBC boundary condition B to initialize.
   * @param bcC DirichletBC boundary condition C to initialize.
   */
  void _initialize(topology::Mesh* mesh,
		   DirichletBC* const bcA,
		   DirichletBC* const bcB,
		   DirichletBC* const bcC) const;

}; // class TestDirichletBCMulti

#endif // pylith_bc_dirichletbcmulti_hh


// End of file 
