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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/faults/TestTractPerturbation.hh
 *
 * @brief C++ TestTractPerturbation object
 *
 * C++ unit testing for TractPerturbation.
 */

#if !defined(pylith_faults_testtractperturbation_hh)
#define pylith_faults_testtractperturbation_hh

#include "pylith/faults/faultsfwd.hh" // USES TractPerturbation
#include "pylith/topology/topologyfwd.hh" // USES Mesh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestTractPerturbation;
  } // faults
} // pylith

/// C++ unit testing for TractPerturbation
class pylith::faults::TestTractPerturbation : public CppUnit::TestFixture
{ // class TestTractPerturbation

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestTractPerturbation );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testLabel );
  CPPUNIT_TEST( testHasParameter );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testCalculate );
  CPPUNIT_TEST( testParameterFields );
  CPPUNIT_TEST( testVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test label().
  void testLabel(void);

  /// Test hasParameter().
  void testHasParameter(void);

  /// Test initialize() with 2-D mesh.
  void testInitialize(void);

  /// Test calculate() with 2-D mesh.
  void testCalculate(void);

  /// Test parameterFields() with 2-D mesh.
  void testParameterFields(void);

  /// Test VertexField() with 2-D mesh.
  void testVertexField(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize TractPerturbation.
   *
   * @param mesh Finite-element mesh of domain.
   * @param faultMesh Finite-element mesh of fault.
   * @param tract Traction perturbation.
   */
  static
  void _initialize(topology::Mesh* mesh,
		   topology::SubMesh* faultMesh,
		   TractPerturbation* tract);

}; // class TestTractPerturbation

#endif // pylith_faults_testtractperturbation_hh


// End of file 
