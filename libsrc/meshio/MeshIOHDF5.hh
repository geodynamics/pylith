// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#if !defined(pylith_meshio_meshiohdf5_hh)
#define pylith_meshio_meshiohdf5_hh

#include <string> // HASA std::string

namespace ALE {
  template<typename T> class Obj;
  class PetscMesh;
} // ALE

namespace pylith {
  namespace meshio {
    class MeshIO;
    class MeshIOHDF5;
    class HDF5; // USES HDF5
  } // pylith
} // meshio

typedef int hid_t;

class pylith::meshio::MeshIOHDF5 : public pylith::meshio::MeshIO
{ // MeshIOHDF5
  
// PUBLIC METHODS -------------------------------------------------------
public :

  /// Constructor
  MeshIOHDF5(void);

  /// Destructor
  ~MeshIOHDF5(void);

  /** Set filename for HDF5 file.
   *
   * @param filename Name of file
   */
  void filename(const char* name);

  /** Get filename of HDF5 file.
   *
   * @returns Name of file
   */
  const char* filename(void) const;

  /** Read mesh from file.
   *
   * @param pMesh Pointer to PETSc mesh object
   */
  void read(ALE::Obj<ALE::PetscMesh>* pMesh);

  /** Write mesh to file.
   *
   * @param mesh PETSc mesh object
   */
  void write(const ALE::Obj<ALE::PetscMesh>& mesh) const;

// PRIVATE METHODS ------------------------------------------------------
private :

  /** Read general mesh information.
   *
   * @param filein HDF5 file
   * @param pMesh Pointer to mesh
   */
  void _readMeshInfo(hid_t& filein,
		     ALE::Obj<ALE::PetscMesh>* pMesh);

  /** Write general mesh information.
   *
   * @param filein HDF5 file
   * @param mesh PETSc mesh
   */
  void _writeMeshInfo(hid_t& fileout,
		     const ALE::Obj<ALE::PetscMesh>& mesh) const;

  /** Read mesh vertices.
   *
   * @param filein HDF5 file
   * @param pCoordinates Pointer to array of vertex coordinates
   * @param pNumVertices Pointer to number of vertices
   * @param pNumDims Pointer to number of dimensions
   */
  void _readVertices(hid_t& filein,
		     double** pCoordinates,
		     int* pNumVertices, 
		     int* pNumDims) const;
  
  /** Write mesh vertices.
   *
   * @param fileout HDF file
   * @param mesh PETSc mesh
   */
  void _writeVertices(hid_t& fileout,
		      const ALE::Obj<ALE::PetscMesh>& mesh) const;
  
  /** Read mesh elements.
   *
   * @param filein Input stream
   * @param pElements Pointer to array of elements (vertices in each element)
   * @param pNumElements Pointer to number of elements
   * @param pNumCorners Pointer to number of corners in each element
   */
  void _readElements(hid_t& filein,
		     int** pElements,
		     int* pNumElements, 
		     int* pNumCorners) const;
  
  /** Write mesh elements.
   *
   * @param fileout HDF file
   * @param mesh PETSc mesh
   */
  void _writeElements(hid_t& fileout,
		      const ALE::Obj<ALE::PetscMesh>& mesh) const;

  /** Read mesh chart.
   *
   * @param filein HDF5 file
   * @param pMesh Pointer to PETSc mesh
   */
  void _readChart(hid_t& filein,
		  ALE::Obj<ALE::PetscMesh>* pMesh) const;

  /** Write mesh chart.
   *
   * @param fileout HDF5 file
   * @param mesh PETSc mesh
   * @param name Name of chart
   */
  void _writeChart(hid_t& fileout,
		   const ALE::Obj<ALE::PetscMesh>& mesh,
		   const char* name) const;

// PRIVATE MEMBERS ------------------------------------------------------
private :

  std::string _filename; ///< Name of file

}; // MeshIOHDF5

#include "MeshIOHDF5.icc" // inline methods

#endif // pylith_meshio_meshiohdf5_hh

// version
// $Id$

// End of file 
