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

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestMeshIO;
    class MeshIO; // USES MeshIO

    class MeshData; // test data
  } // meshio

  namespace topology {
    class Mesh; // USES Mesh
  } // topology
} // pylith

// MeshIO ---------------------------------------------------------------
class pylith::meshio::TestMeshIO : public CppUnit::TestFixture
{ // class TestMeshIO

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Get simple mesh for testing I/O.
   *
   * @param data Mesh data
   *
   * @returns PyLith mesh
   */
  topology::Mesh* _createMesh(const MeshData& data);

  /** Check values in mesh against data.
   *
   * @param mesh PyLith mesh
   * @param data Mesh data
   */
  void _checkVals(const topology::Mesh& mesh,
		  const MeshData& data);

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

}; // class TestMeshIO

#endif // pylith_meshio_testmeshio_hh


// End of file 
