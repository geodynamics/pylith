// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/feassemble/TestElasticityExplicitTri3.hh
 *
 * @brief C++ TestElasticityExplicitTri3 object
 *
 * C++ unit testing for ElasticityExplicitTri3.
 */

#if !defined(pylith_feassemble_testelasticityexplicittri3_hh)
#define pylith_feassemble_testelasticityexplicittri3_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/feassemble/feassemblefwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // USES Mesh, SolutionFields
#include "pylith/materials/materialsfwd.hh" // USES ElasticMaterial

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityExplicitTri3;
    class ElasticityExplicitData;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityExplicitTri3
class pylith::feassemble::TestElasticityExplicitTri3 : public CppUnit::TestFixture
{ // class TestElasticityExplicitTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityExplicitTri3 );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testMaterial );
  CPPUNIT_TEST( testNeedNewJacobian );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
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

  /// Test initialize().
  void testInitialize(void);

  /// Test integrateResidual().
  void testIntegrateResidual(void);

  /// Test integrateJacobian().
  void testIntegrateJacobian(void);

  /// Test updateStateVars().
  void testUpdateStateVars(void);

  /// Test StableTimeStep().
  void testStableTimeStep(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  ElasticityExplicitData* _data; ///< Data for testing.
  materials::ElasticMaterial* _material; ///< Elastic material.
  Quadrature* _quadrature; ///< Quadrature information.
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
		   ElasticityExplicitTri3* const integrator,
		   topology::SolutionFields* const fields);

}; // class TestElasticityExplicitTri3

#endif // pylith_feassemble_testelasticityexplicittri3_hh


// End of file 
