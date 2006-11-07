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

#include "MeshIO.hh"

namespace pylith {
  namespace meshio {
    class MeshIOAscii;
  } // meshio
} // pylith

class pylith::meshio::MeshIOAscii : public MeshIO
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

// PROTECTED METHODS ----------------------------------------------------
protected :

  /// Write mesh
  void _write(void) const;

  /// Read mesh
  void _read(void);

// PRIVATE METHODS ------------------------------------------------------
private :

  /** Read mesh vertices.
   *
   * @param filein Input stream
   * @param pCoordinates Pointer to array of vertex coordinates
   * @param pNumVertices Pointer to number of vertices
   * @param pSpaceDim Pointer to dimension of coordinates vector space
   */
  void _readVertices(std::istream& filein,
		     double** pCoordinates,
		     int* pNumVertices,
		     int* pSpaceDim) const;
  
  /** Write mesh vertices.
   *
   * @param fileout Output stream
   * @param coordinates Array of vertex coordinates
   * @param numVertices Number of vertices
   * @param spaceDim Dimension of coordinates vector space
   */
  void _writeVertices(std::ostream& fileout,
		      const double* coordinates,
		      const int numVertices,
		      const int spaceDim) const;
  
  /** Read mesh cells.
   *
   * @param filein Input stream
   * @param pCells Pointer to array of indices of cell vertices
   * @param pNumCells Pointer to number of cells
   * @param pNumCorners Pointer to number of corners
   */
  void _readCells(std::istream& filein,
		  int** pCells,
		  int* pNumCells,
		  int* pNumCorners) const;
  
  /** Write mesh cells.
   *
   * @param fileout Output stream
   * @param cells Array of indices of cell vertices
   * @param numCells Number of cells
   * @param numCorners Number of corners
   */
  void _writeCells(std::ostream& fileout,
		   const int* cells,
		   const int numCells,
		   const int numCorners) const;

  // PRIVATE MEMBERS ----------------------------------------------------
private :

  std::string _filename; ///< Name of file

}; // MeshIOAscii

#include "MeshIOAscii.icc" // inline methods

#endif // pylith_meshio_meshioascii_hh

// End of file 
