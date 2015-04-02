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
 * @file unittests/libtests/meshio/TestMeshIO.hh
 *
 * @brief C++ TestMeshIO object
 *
 * C++ unit testing for MeshIO.
 */

#if !defined(pylith_meshio_testmeshio_hh)
#define pylith_meshio_testmeshio_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "pylith/meshio/meshiofwd.hh" // USES MeshIO

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestMeshIO;

    class MeshData; // test data
  } // meshio
} // pylith

// MeshIO ---------------------------------------------------------------
class pylith::meshio::TestMeshIO : public CppUnit::TestFixture
{ // class TestMeshIO

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Get simple mesh for testing I/O.
   *
   * @param data Mesh data
   *
   * @returns PyLith mesh
   */
  void _createMesh(const MeshData& data);

  /** Check values in mesh against data.
   *
   * @param data Mesh data
   */
  void _checkVals(const MeshData& data);

  /** Test debug().
   *
   * @param iohandler MeshIO object.
   */
  void _testDebug(MeshIO& iohandler);

  /** Test interpolate().
   *
   * @param iohandler MeshIO object.
   */
  void _testInterpolate(MeshIO& iohandler);

// PROTECTED MEMEBERS ////////////////////////////////////////////////////
protected :

  topology::Mesh* _mesh; ///< Finite-element mesh.

}; // class TestMeshIO

#endif // pylith_meshio_testmeshio_hh


// End of file 
