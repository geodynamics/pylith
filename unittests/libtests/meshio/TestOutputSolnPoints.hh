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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/meshio/TestOutputSolnPoints.hh
 *
 * @brief C++ TestOutputSolnPoints object
 *
 * C++ unit testing for OutputSolnPoints.
 */

#if !defined(pylith_meshio_testoutputsolnpoints_hh)
#define pylith_meshio_testoutputsolnpoints_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // USES Field

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestOutputSolnPoints;

    class OutputSolnPointsData;
  } // meshio
} // pylith

/// C++ unit testing for OutputSolnPoints
class pylith::meshio::TestOutputSolnPoints : public CppUnit::TestFixture
{ // class TestOutputSolnPoints

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE( TestOutputSolnPoints );
    
    CPPUNIT_TEST( testConstructor );
    
    CPPUNIT_TEST( testSetupInterpolatorTri3 );
    CPPUNIT_TEST( testInterpolateTri3 );

    CPPUNIT_TEST( testSetupInterpolatorQuad4 );
    CPPUNIT_TEST( testInterpolateQuad4 );

    CPPUNIT_TEST( testSetupInterpolatorTet4 );
    CPPUNIT_TEST( testInterpolateTet4 );

    CPPUNIT_TEST( testSetupInterpolatorHex8 );
    CPPUNIT_TEST( testInterpolateHex8 );

    CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test setupInterpolator for tri3 mesh.
  void testSetupInterpolatorTri3(void);

  /// Test interpolation for tri3 mesh.
  void testInterpolateTri3(void);

  /// Test setupInterpolator for quad4 mesh.
  void testSetupInterpolatorQuad4(void);

  /// Test interpolation for quad4 mesh.
  void testInterpolateQuad4(void);

  /// Test setupInterpolator for tet4 mesh.
  void testSetupInterpolatorTet4(void);

  /// Test interpolation for tet4 mesh.
  void testInterpolateTet4(void);

  /// Test setupInterpolator for hex8 mesh.
  void testSetupInterpolatorHex8(void);

  /// Test interpolation for hex8 mesh.
  void testInterpolateHex8(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Test setupInterpolator.
   *
   * @param data Test data.
   */
  void _testSetupInterpolator(const OutputSolnPointsData& data);

  /** Test interpolation.
   *
   * @param data Test data.
   */
  void _testInterpolate(const OutputSolnPointsData& data);

  /** Compute values of field at vertices in mesh.
   *
   * @param field Field to hold values.
   * @param data Test data.
   */
  void _calcField(pylith::topology::Field* field,
		  const OutputSolnPointsData& data);

}; // class TestOutputSolnPoints

#endif // pylith_meshio_testoutputsolnpoints_hh

// End of file 
