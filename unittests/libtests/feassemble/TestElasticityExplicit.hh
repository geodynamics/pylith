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
 * @file unittests/libtests/feassemble/TestElasticityExplicit.hh
 *
 * @brief C++ TestElasticityExplicit object
 *
 * C++ unit testing for ElasticityExplicit.
 */

#if !defined(pylith_feassemble_testelasticityexplicit_hh)
#define pylith_feassemble_testelasticityexplicit_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityExplicit;

    class ElasticityExplicit; // USES ElasticityExplicit
    class IntegratorData; // HOLDSA IntegratorData
    class Quadrature; // HOLDSA Quadrature
  } // feassemble

  namespace materials {
    class ElasticMaterial; // HOLDSA ElasticMaterial
  } // materials

  namespace topology {
    class FieldsManager; // USES FieldsManager
  } // topology
} // pylith

/// C++ unit testing for ElasticityExplicit
class pylith::feassemble::TestElasticityExplicit : public CppUnit::TestFixture
{ // class TestElasticityExplicit

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicit );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testStableTimeStep );
  CPPUNIT_TEST( testMaterial );
  CPPUNIT_TEST( testNeedNewJacobian );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test timeStep().
  void testTimeStep(void);

  /// Test StableTimeStep().
  void testStableTimeStep(void);

  /// Test material().
  void testMaterial(void);

  /// Test needNewJacobian().
  void testNeedNewJacobian(void);

  /// Test updateState().
  void testUpdateState(void);

  /// Test integrateResidual().
  void testIntegrateResidual(void);

  /// Test integrateJacobian().
  void testIntegrateJacobian(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  IntegratorData* _data; ///< Data for testing.
  materials::ElasticMaterial* _material; ///< Elastic material.
  Quadrature* _quadrature; ///< Quadrature information.

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize elasticity integrator.
   *
   * @param mesh PETSc mesh to initialize.
   * @param integrator ElasticityIntegrator to initialize.
   * @param fields Solution fields.
   */
  void _initialize(ALE::Obj<ALE::Mesh>* mesh,
		   ElasticityExplicit* const integrator,
		   topology::FieldsManager* const fields);

}; // class TestElasticityExplicit

#endif // pylith_feassemble_testelasticityexplicit_hh


// End of file 
