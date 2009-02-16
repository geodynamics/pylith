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
 * @file unittests/libtests/bc/TestNeumann.hh
 *
 * @brief C++ TestNeumann object.
 *
 * C++ unit testing for Neumann.
 */

#if !defined(pylith_bc_testneumann_hh)
#define pylith_bc_testneumann_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/bc/bcfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/feassemble/feassemblefwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestNeumann;
    class NeumannData; // HOLDSA NeumannData
  } // bc
} // pylith

/// C++ unit testing for Neumann.
class pylith::bc::TestNeumann : public CppUnit::TestFixture
{ // class TestNeumann

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestNeumann );
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

  /// Test integrateResidual().
  void testIntegrateResidual(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  NeumannData* _data; ///< Data for testing
  feassemble::Quadrature<topology::SubMesh>* _quadrature; ///< Used in testing.

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize Neumann boundary condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param bc Neumann boundary condition to initialize.
   * @param fields Solution fields.
   */
  void _initialize(topology::Mesh* mesh,
		   Neumann* const bc,
		   topology::SolutionFields* fields) const;

}; // class TestNeumann

#endif // pylith_bc_neumann_hh


// End of file 
