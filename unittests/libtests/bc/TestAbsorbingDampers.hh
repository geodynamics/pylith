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
 * @file unittests/libtests/bc/TestAbsorbingDampers.hh
 *
 * @brief C++ TestAbsorbingDampers object.
 *
 * C++ unit testing for AbsorbingDampers.
 */

#if !defined(pylith_bc_testabsorbingdampers_hh)
#define pylith_bc_testabsorbingdampers_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestAbsorbingDampers;

    class AbsorbingDampers;
    class AbsorbingDampersData;
  } // bc

  namespace feassemble {
    class Quadrature; // HOLDSA Quadrature
  } // feassemble

  namespace topology {
    class FieldsManager; // USES FieldsManager
  } // topology
} // pylith

/// C++ unit testing for AbsorbingDampers.
class pylith::bc::TestAbsorbingDampers : public CppUnit::TestFixture
{ // class TestAbsorbingDampers

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestAbsorbingDampers );
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

  /// Test integrateJacobian().
  void testIntegrateJacobian(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  AbsorbingDampersData* _data; ///< Data for testing
  feassemble::Quadrature* _quadrature; ///< Data used in testing

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize AbsorbingDampers boundary condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param bc AbsorbingDampers boundary condition to initialize.
   * @param cs Mesh coordinate system.
   * @param fields Solution fields.
   */
  void _initialize(ALE::Obj<Mesh>* mesh,
		   AbsorbingDampers* const bc,
		   topology::FieldsManager* fields) const;

}; // class TestAbsorbingDampers

#endif // pylith_bc_absorbingdampers_hh


// End of file 
