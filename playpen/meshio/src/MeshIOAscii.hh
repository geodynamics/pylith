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

#if !defined(meshioascii_hh)
#define meshioascii_hh

#include <iosfwd> // USES std::istream, std::ostream
#include <string> // HASA std::string

namespace ALE {
  template<typename T> class Obj;
  class PetscMesh;
}

class MeshIOAscii
{ // MeshIOAscii
  
// PUBLIC METHODS -------------------------------------------------------
public :

  /// Constructor
  MeshIOAscii(void);

  /// Destructor
  ~MeshIOAscii(void);

  /** Set filename for ASCII file.
   *
   * @param filename Name of file
   */
  void filename(const char* name);

  /** Get filename of ASCII file.
   *
   * @returns Name of file
   */
  const char* filename(void) const;

  /** Read mesh from file.
   *
   * @param pMesh Pointer to PETSc mesh object
   */
  void read(ALE::Obj<ALE::PetscMesh>* pMesh) const;

  /** Write mesh to file.
   *
   * @param mesh PETSc mesh object
   */
  void write(const ALE::Obj<ALE::PetscMesh>& mesh) const;

// PRIVATE METHODS ------------------------------------------------------
private :

  /** Read mesh vertices.
   *
   * @param filein Input stream
   * @param pCoordinates Pointer to array of vertex coordinates
   * @param pNumVertices Pointer to number of vertices
   * @param pNumDims Pointer to number of dimensions
   */
  void _readVertices(std::istream& filein,
		     double** pCoordinates,
		     int* pNumVertices, 
		     int* pNumDims) const;
  
  /** Write mesh vertices.
   *
   * @param fileout Output stream
   * @param mesh PETSc mesh
   */
  void _writeVertices(std::ostream& fileout,
		      const ALE::Obj<ALE::PetscMesh>& mesh) const;
  
  /** Read mesh elements.
   *
   * @param filein Input stream
   * @param pElements Pointer to array of elements (vertices in each element)
   * @param pNumElements Pointer to number of elements
   * @param pNumCorners Pointer to number of corners in each element
   */
  void _readElements(std::istream& filein,
		     int** pElements,
		     int* pNumElements, 
		     int* pNumCorners) const;
  
  /** Write mesh elements.
   *
   * @param fileout Output stream
   * @param mesh PETSc mesh
   */
  void _writeElements(std::ostream& fileout,
		      const ALE::Obj<ALE::PetscMesh>& mesh) const;

  /** Read mesh chart.
   *
   * @param filein Output stream
   * @param pMesh Pointer to PETSc mesh
   */
  void _readChart(std::istream& filein,
		  ALE::Obj<ALE::PetscMesh>* pMesh) const;

  /** Write mesh chart.
   *
   * @param fileout Output stream
   * @param mesh PETSc mesh
   * @param name Name of chart
   */
  void _writeChart(std::ostream& fileout,
		   const ALE::Obj<ALE::PetscMesh>& mesh,
		   const char* name) const;

// PRIVATE MEMBERS ------------------------------------------------------
private :

  std::string _filename; ///< Name of file

}; // MeshIOAscii

#include "MeshIOAscii.icc" // inline methods

#endif // meshioascii_hh

// version
// $Id$

// End of file 
