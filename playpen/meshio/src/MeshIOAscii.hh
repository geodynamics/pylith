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

#if !defined(pylith_meshio_meshioascii_hh)
#define pylith_meshio_meshioascii_hh

#include <iosfwd> // USES std::istream, std::ostream
#include <string> // HASA std::string
#include <MeshIO.hh>

using ALE::Obj;

namespace pylith {
  namespace meshIO {
    class MeshIOAscii;
  } // meshio
} // pylith

class pylith::meshIO::MeshIOAscii : public pylith::meshIO::MeshIO
{ // MeshIOAscii
// PUBLIC TYPEDEFS -------------------------------------------------------
public :
  typedef ALE::Mesh                Mesh;
  typedef ALE::Mesh::sieve_type    sieve_type;
  typedef ALE::Mesh::topology_type topology_type;
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
  void read(Obj<Mesh>& mesh, const bool interpolate = false);

  /** Write mesh to file.
   *
   * @param mesh PETSc mesh object
   */
  void write(const Obj<Mesh>& mesh) const;

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
		      const Obj<Mesh>& mesh) const;
  
  /** Read mesh cells.
   *
   * @param filein Input stream
   * @param pCells Pointer to array of cells (vertices in each element)
   * @param pNumCells Pointer to number of cells
   * @param pNumCorners Pointer to number of corners in each element
   */
  void _readCells(std::istream& filein,
		     int** pCells,
		     int* pNumCells, 
		     int* pNumCorners) const;
  
  /** Write mesh cells.
   *
   * @param fileout Output stream
   * @param mesh PETSc mesh
   */
  void _writeCells(std::ostream& fileout,
		      const Obj<Mesh>& mesh) const;

  /** Read mesh chart.
   *
   * @param filein Output stream
   * @param pMesh Pointer to PETSc mesh
   */
  void _readChart(std::istream& filein,
		  const Obj<Mesh>& pMesh) const;

  /** Write mesh chart.
   *
   * @param fileout Output stream
   * @param mesh PETSc mesh
   * @param name Name of chart
   */
  void _writeChart(std::ostream& fileout,
		   const Obj<Mesh>& mesh,
		   const char* name) const;

// PRIVATE MEMBERS ------------------------------------------------------
private :

  std::string _filename; ///< Name of file

}; // MeshIOAscii

#include "MeshIOAscii.icc" // inline methods

#endif // pylith_meshio_meshioascii_hh

// version
// $Id$

// End of file 
