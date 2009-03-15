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
 * @file unittests/libtests/feassemble/TestElasticityImplicit.hh
 *
 * @brief C++ TestElasticityImplicit object
 *
 * C++ unit testing for ElasticityImplicit.
 */

#if !defined(pylith_feassemble_testelasticityimplicit_hh)
#define pylith_feassemble_testelasticityimplicit_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityImplicit;

    class ElasticityImplicit; // USES ElasticityImplicit
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

namespace spatialdata {
  namespace spatialdb {
    class GravityField; // HOLDSA GravityField
  } // spatialdb
} // spatialdata

/// C++ unit testing for ElasticityImplicit
class pylith::feassemble::TestElasticityImplicit : public CppUnit::TestFixture
{ // class TestElasticityImplicit

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicit );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testStableTimeStep );
  CPPUNIT_TEST( testMaterial );
  CPPUNIT_TEST( testNeedNewJacobian );
  CPPUNIT_TEST( testUseSolnIncr );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

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

  /// Test useSolnIncr().
  void testUseSolnIncr(void);

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
  spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field.

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize elasticity integrator.
   *
   * @param mesh PETSc mesh to initialize.
   * @param integrator ElasticityIntegrator to initialize.
   * @param fields Solution fields.
   */
  void _initialize(ALE::Obj<Mesh>* mesh,
		   ElasticityImplicit* const integrator,
		   topology::FieldsManager* const fields);

}; // class TestElasticityImplicit

#endif // pylith_feassemble_testelasticityimplicit_hh


// End of file 
