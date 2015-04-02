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
 * @file unittests/libtests/materials/TestMaterial.hh
 *
 * @brief C++ TestMaterial object
 *
 * C++ unit testing for Material.
 */

#if !defined(pylith_materials_testmaterial_hh)
#define pylith_materials_testmaterial_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/materials/materialsfwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class TestMaterial;
    class MaterialData;
  } // materials
} // pylith

/// C++ unit testing for Material
class pylith::materials::TestMaterial : public CppUnit::TestFixture
{ // class TestMaterial

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMaterial );

  CPPUNIT_TEST( testID );
  CPPUNIT_TEST( testLabel );
  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testDBProperties );
  CPPUNIT_TEST( testDBStateVars );
  CPPUNIT_TEST( testNormalizer );
  CPPUNIT_TEST( testNeedNewJacobian );
  CPPUNIT_TEST( testIsJacobianSymmetric );
  CPPUNIT_TEST( testInitialize );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test id().
  void testID(void);

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

  /// Test needNewJacobian()
  void testNeedNewJacobian(void);

  /// Test isJacobianSymmetric()
  void testIsJacobianSymmetric(void);

  /// Test initialize()
  void testInitialize(void);

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  // Methods used in testing children of this class.

  /// Setup testing data.
  virtual
  void setUp(void);

  /// Tear down testing data.
  virtual
  void tearDown(void);

  /// Test dimension().
  void testDimension();

  /// Test tensorSize().
  void testTensorSize();

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

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  Material* _material; ///< Object for testing
  MaterialData* _data; ///< Data for testing

}; // class TestMaterial

#endif // pylith_materials_testmaterial_hh

// End of file 
