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

#include <Mesh.hh>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestMeshIO;
    class MeshData;
  } // meshio
} // pylith

/// C++ unit testing for TestMeshIO
class pylith::meshio::TestMeshIO : public CppUnit::TestFixture
{ // class TestMeshIO
  // PUBLIC TYPEDEFS ////////////////////////////////////////////////////
public :

  typedef ALE::Field::Mesh Mesh;
  typedef Mesh::sieve_type sieve_type;

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Get simple mesh for testing I/O.
   *
   * @param data Mesh data
   *
   * @returns PETSc mesh
   */
  ALE::Obj<Mesh>* createMesh(const MeshData& data);

  /** Check values in mesh against data.
   *
   * @param mesh PETSc mesh
   * @param data Mesh data
   */
  void checkVals(const ALE::Obj<Mesh>& mesh,
		 const MeshData& data);

}; // class TestMeshIO

#endif // pylith_meshio_testmeshio_hh

// End of file 
