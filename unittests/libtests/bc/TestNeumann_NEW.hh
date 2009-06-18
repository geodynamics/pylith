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
 * @file unittests/libtests/bc/TestNeumann_NEW.hh
 *
 * @brief C++ TestNeumann_NEW object.
 *
 * C++ unit testing for Neumann_NEW.
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
    class TestNeumann_NEW;
    class NeumannData; // HOLDSA NeumannData
  } // bc
} // pylith

/// C++ unit testing for Neumann_NEW.
class pylith::bc::TestNeumann_NEW : public CppUnit::TestFixture
{ // class TestNeumann_NEW

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestNeumann_NEW );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( test_getLabel );
  CPPUNIT_TEST( test_queryDB );
  CPPUNIT_TEST( test_queryDatabases );
  CPPUNIT_TEST( test_paramsLocalToGlobal );
  CPPUNIT_TEST( test_calculateValueInitial );
  CPPUNIT_TEST( test_calculateValueRate );
  CPPUNIT_TEST( test_calculateValueChange );
  CPPUNIT_TEST( test_calculateValueChangeTH );
  CPPUNIT_TEST( test_calculateValueAll );

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

  /// Test _getLabel().
  void test_getLabel(void);

  /// Test _queryDB().
  void test_queryDB(void);

  /// Test _queryDatabases().
  void test_queryDatabases(void);

  /// Test _paramsLocalToGlobal().
  void test_paramsLocalToGlobal(void);

  /// Test _calculateValue() with initial value.
  void test_calculateValueInitial(void);

  /// Test _calculateValue() with rate.
  void test_calculateValueRate(void);

  /// Test _calculateValue() with temporal change.
  void test_calculateValueChange(void);

  /// Test _calculateValue() with temporal change w/time history.
  void test_calculateValueChangeTH(void);

  /// Test _calculateValue() with initial, rate, and temporal change
  /// w/time history.
  void test_calculateValueAll(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  NeumannData* _data; ///< Data for testing
  feassemble::Quadrature<topology::SubMesh>* _quadrature; ///< Used in testing.

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Do minimal initialization of Neumann_NEW boundary condition.
   *
   * @param mesh Finite-element mesh to initialize
   * @param bc Neumann_NEW boundary condition to initialize.
   * @param useScales Use scales provided by local constants.
   */
  void _preinitialize(topology::Mesh* mesh,
		      Neumann_NEW* const bc,
		      const bool useScales =false) const;

  /** Initialize Neumann_NEW boundary condition.
   *
   * @param mesh Finite-element mesh to initialize
   * @param bc Neumann_NEW boundary condition to initialize.
   * @param fields Solution fields.
   */
  void _initialize(topology::Mesh* mesh,
		   Neumann_NEW* const bc,
		   topology::SolutionFields* fields) const;

}; // class TestNeumann_NEW

#endif // pylith_bc_neumann_hh


// End of file 
