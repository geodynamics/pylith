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
 * @file unittests/libtests/bc/TestAbsorbingDampersQuad4.hh
 *
 * @brief C++ TestAbsorbingDampers object.
 *
 * C++ unit testing for AbsorbingDampers for mesh with 2-D quad cells.
 */

#if !defined(pylith_bc_testabsorbingdampersquad4_hh)
#define pylith_bc_testabsorbingdampersquad4_hh

#include "TestAbsorbingDampers.hh" // ISA TestAbsorbingDampers

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestAbsorbingDampersQuad4;
  } // bc
} // pylith

/// C++ unit testing for AbsorbingDampers for mesh with 2-D quad cells.
class pylith::bc::TestAbsorbingDampersQuad4 : public TestAbsorbingDampers
{ // class TestAbsorbingDampers

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestAbsorbingDampersQuad4, TestAbsorbingDampers );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestAbsorbingDampersQuad4

#endif // pylith_bc_absorbingdampersquad4_hh


// End of file 
