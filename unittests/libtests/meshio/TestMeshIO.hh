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

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestMeshIO;
    class MeshIO;

    class MeshData;
  } // meshio
} // pylith

/// C++ unit testing for TestMeshIO
class pylith::meshio::TestMeshIO : public CppUnit::TestFixture
{ // class TestMeshIO

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Get simple mesh for testing I/O.
   *
   * @param data Mesh data
   *
   * @returns PETSc mesh
   */
  ALE::Obj<ALE::Mesh>* _createMesh(const MeshData& data);

  /** Check values in mesh against data.
   *
   * @param mesh PETSc mesh
   * @param data Mesh data
   */
  void _checkVals(const ALE::Obj<ALE::Mesh>& mesh,
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
