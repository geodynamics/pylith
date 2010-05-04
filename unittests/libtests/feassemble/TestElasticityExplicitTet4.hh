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
 * @file unittests/libtests/feassemble/TestElasticityExplicitTet4.hh
 *
 * @brief C++ TestElasticityExplicitTet4 object
 *
 * C++ unit testing for ElasticityExplicitTet4.
 */

#if !defined(pylith_feassemble_testelasticityexplicittet4_hh)
#define pylith_feassemble_testelasticityexplicittet4_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/feassemble/feassemblefwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // USES Mesh, SolutionFields
#include "pylith/materials/materialsfwd.hh" // USES ElasticMaterial

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityExplicitTet4;
    class ElasticityExplicitData;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityExplicitTet4
class pylith::feassemble::TestElasticityExplicitTet4 : public CppUnit::TestFixture
{ // class TestElasticityExplicitTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitTet4 );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testMaterial );
  CPPUNIT_TEST( testNeedNewJacobian );
  CPPUNIT_TEST( testUseSolnIncr );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateResidualLumped );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testStableTimeStep );

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

  /// Test material().
  void testMaterial(void);

  /// Test needNewJacobian().
  void testNeedNewJacobian(void);

  /// Test useSolnIncr().
  void testUseSolnIncr(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test integrateResidual().
  void testIntegrateResidual(void);

  /// Test integrateResidualLumped().
  void testIntegrateResidualLumped(void);

  /// Test integrateJacobian().
  void testIntegrateJacobian(void);

  /// Test integrateJacobianLumped().
  void testIntegrateJacobianLumped(void);

  /// Test updateStateVars().
  void testUpdateStateVars(void);

  /// Test StableTimeStep().
  void testStableTimeStep(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  ElasticityExplicitData* _data; ///< Data for testing.
  materials::ElasticMaterial* _material; ///< Elastic material.
  Quadrature<topology::Mesh>* _quadrature; ///< Quadrature information.
  spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field.

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize elasticity integrator.
   *
   * @param mesh Finite-element mesh to initialize.
   * @param integrator ElasticityIntegrator to initialize.
   * @param fields Solution fields.
   */
  void _initialize(topology::Mesh* mesh,
		   ElasticityExplicitTet4* const integrator,
		   topology::SolutionFields* const fields);

}; // class TestElasticityExplicitTet4

#endif // pylith_feassemble_testelasticityexplicittet4_hh


// End of file 
