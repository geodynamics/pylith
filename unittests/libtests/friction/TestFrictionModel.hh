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
 * @file unittests/libtests/friction/TestFrictionModel.hh
 *
 * @brief C++ TestFrictionModel object
 *
 * C++ unit testing for FrictionModel.
 */

#if !defined(pylith_friction_testfrictionmodel_hh)
#define pylith_friction_testfrictionmodel_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/friction/frictionfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/faults/faultsfwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace friction {
    class TestFrictionModel;
    class FrictionModelData;
    class StaticFrictionData;
  } // friction
} // pylith

/// C++ unit testing for FrictionModel
class pylith::friction::TestFrictionModel : public CppUnit::TestFixture
{ // class TestFrictionModel

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFrictionModel );

  CPPUNIT_TEST( testLabel );
  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testDBProperties );
  CPPUNIT_TEST( testDBStateVars );
  CPPUNIT_TEST( testNormalizer );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testGetField );
  CPPUNIT_TEST( testRetrievePropsStateVars );
  CPPUNIT_TEST( testCalcFriction );
  CPPUNIT_TEST( testCalcFrictionDeriv );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test label()
  void testLabel(void);

  /// Test timeStep()
  void testTimeStep(void);

  /// Test dbProperties()
  void testDBProperties(void);

  /// Test dbStateVars().
  void testDBStateVars(void);

  /// Test normalizer().
  void testNormalizer(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test getField().
  void testGetField(void);

  /// Test retrievePropsStateVars().
  void testRetrievePropsStateVars(void);

  /// Test calcFriction()
  void testCalcFriction(void);

  /// Test calcFrictionDeriv()
  void testCalcFrictionDeriv(void);

  /// Test updateStateVars().
  void testUpdateStateVars(void);

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  // Methods used in testing children of this class.

  /// Setup testing data.
  virtual
  void setUp(void);

  /// Tear down testing data.
  virtual
  void tearDown(void);

  /// Test _dbToProperties().
  void testDBToProperties(void);

  /// Test _nondimProperties().
  void testNonDimProperties(void);

  /// Test _dimProperties().
  void testDimProperties(void);

  /// Test _dbToStateVars().
  void testDBToStateVars(void);

  /// Test _nondimStateVars().
  void testNonDimStateVars(void);

  /// Test _dimStateVars().
  void testDimStateVars(void);

  /// Test _calcFriction().
  void test_calcFriction(void);

  /// Test _calcFrictionDeriv().
  void test_calcFrictionDeriv(void);

  /// Test _updateStateVars().
  void test_updateStateVars(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /// Setup nondimensionalization.
  void setupNormalizer(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  FrictionModel* _friction; ///< Test subject.
  FrictionModelData* _data; ///< Data for tests.

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /** Setup mesh and material.
   *
   * @param mesh Finite-element mesh.
   * @param fault Fault with friction.
   * @param friction Friction model.
   * @param data Data with properties for friction model.
   */
  void _initialize(topology::Mesh* mesh,
                   faults::FaultCohesiveDyn* fault,
                   StaticFriction* friction,
                   const StaticFrictionData* data);

}; // class TestFrictionModel

#endif // pylith_friction_testfrictionmodel_hh

// End of file
