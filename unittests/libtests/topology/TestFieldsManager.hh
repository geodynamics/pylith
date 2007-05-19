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
 * @file unittests/libtests/topology/TestFieldsManager.hh
 *
 * @brief C++ TestFieldsManager object.
 *
 * C++ unit testing for FieldsManager.
 */

#if !defined(pylith_topology_testfieldsmanager_hh)
#define pylith_topology_testfieldsmanager_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestFieldsManager;

    class FieldsManager;
  } // topology
} // pylith

/// C++ unit testing for FieldsManager.
class pylith::topology::TestFieldsManager : public CppUnit::TestFixture
{ // class TestFieldsManager

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFieldsManager );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testAddReal );
  CPPUNIT_TEST( testGetReal );
  CPPUNIT_TEST( testDelReal );
  CPPUNIT_TEST( testSetFiberDimension );
  CPPUNIT_TEST( testAllocate );
  CPPUNIT_TEST( testCopyLayout );
  CPPUNIT_TEST( testCopyLayoutFromField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test addReal().
  void testAddReal(void);

  /// Test getReal().
  void testGetReal(void);

  /// Test delReal().
  void testDelReal(void);

  /// Test setFiberDimension().
  void testSetFiberDimension(void);

  /// Test allocate().
  void testAllocate(void);

  /// Test copyLayout().
  void testCopyLayout(void);

  /// Test copyLayoutFromField().
  void testCopyLayoutFromField(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize FieldsManager boundary condition.
   *
   * @param mesh PETSc mesh to initialize
   */
  void _initialize(ALE::Obj<ALE::Mesh>* mesh) const;

}; // class TestFieldsManager

#endif // pylith_topology_fieldsmanager_hh


// End of file 
