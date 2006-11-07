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
 * @file unittests/libtests/feassemble/TestIntegratorInertia.hh
 *
 * @brief C++ TestIntegratorInertia object
 *
 * C++ unit testing for IntegratorInertia.
 */

#if !defined(pylith_feassemble_testintegratorinertia_hh)
#define pylith_feassemble_testintegratorinertia_hh

#include "TestIntegrator.hh"

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestIntegratorInertia;

    class IntegratorInertia; // USES IntegratorInertia
  } // feassemble
} // pylith

/// C++ unit testing for IntegratorInertia
class pylith::feassemble::TestIntegratorInertia : public TestIntegrator
{ // class TestIntegratorInertia1D

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Test integrateLumped()
   *
   * @param integrator Pointer to integrator
   * @param data Data for testing integrator
   */
  void _testIntegrateLumped(IntegratorInertia* integrator,
			    const IntegratorData& data) const;

}; // class TestIntegratorInertia

#endif // pylith_feassemble_testintegratorinertia_hh

// End of file 
