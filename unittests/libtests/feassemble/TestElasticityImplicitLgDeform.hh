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
 * @file unittests/libtests/feassemble/TestElasticityImplicitLgDeform.hh
 *
 * @brief C++ TestElasticityImplicitLgDeform object
 *
 * C++ unit testing for ElasticityImplicitLgDeform.
 */

#if !defined(pylith_feassemble_testelasticityimplicitlgdeform_hh)
#define pylith_feassemble_testelasticityimplicitlgdeform_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/feassemble/feassemblefwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // USES Mesh, SolutionFields
#include "pylith/materials/materialsfwd.hh" // USES ElasticMaterial

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestElasticityImplicitLgDeform;
    class IntegratorData;
  } // feassemble
} // pylith

/// C++ unit testing for ElasticityImplicitLgDeform
class pylith::feassemble::TestElasticityImplicitLgDeform : public CppUnit::TestFixture
{ // class TestElasticityImplicitLgDeform

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticityImplicitLgDeform );

  CPPUNIT_TEST( testConstructor );

  // Testing of initialize(), integrateResidual(),
  // integrateJacobian(), and updateStateVars() handled by derived
  // classes.

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

  /// Test updateStateVars().
  void testUpdateStateVars(void);

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
   * @param mesh Finite-element mesh to initialize.
   * @param integrator ElasticityIntegrator to initialize.
   * @param fields Solution fields.
   */
  void _initialize(topology::Mesh* mesh,
		   ElasticityImplicitLgDeform* const integrator,
		   topology::SolutionFields* const fields);

}; // class TestElasticityImplicitLgDeform

#endif // pylith_feassemble_testelasticityimplicitlgdeform_hh


// End of file 
