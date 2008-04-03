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
 * @file unittests/libtests/bc/TestDirichletPointsMulti.hh
 *
 * @brief C++ TestDirichletPointsMulti object.
 *
 * C++ unit testing for DirichletPointsMulti.
 */

#if !defined(pylith_bc_testdirichletpointsmulti_hh)
#define pylith_bc_testdirichletpointsmulti_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletPointsMulti;

    class DirichletPoints;
    class DirichletDataMulti;
  } // bc
} // pylith

/// C++ unit testing for DirichletPointsMulti.
class pylith::bc::TestDirichletPointsMulti : public CppUnit::TestFixture
{ // class TestDirichletPointsMulti

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

  /** Initialize DirichletPointsMulti boundary condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param bcA DirichletPoints boundary condition A to initialize.
   * @param bcB DirichletPoints boundary condition B to initialize.
   * @param bcC DirichletPoints boundary condition C to initialize.
   */
  void _initialize(ALE::Obj<Mesh>* mesh,
		   DirichletPoints* const bcA,
		   DirichletPoints* const bcB,
		   DirichletPoints* const bcC) const;

}; // class TestDirichletPointsMulti

#endif // pylith_bc_dirichletpointsmulti_hh


// End of file 
