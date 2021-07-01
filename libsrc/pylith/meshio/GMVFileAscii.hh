// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/GMVFileAscii.hh
 *
 * @brief C++ object for input/output of LaGriT ASCII GMV files.
 */

#if !defined(pylith_meshio_gmvfileascii_hh)
#define pylith_meshio_gmvfileascii_hh

// Include directives ---------------------------------------------------
#include "GMVFile.hh" // ISA GMVFile

#include "pylith/utils/arrayfwd.hh" // USES int_array, scalar_array
#include <iosfwd>

// GMVFileAscii ---------------------------------------------------------
/// C++ object for input/output of LaGriT ASCII GMV files.
class pylith::meshio::GMVFileAscii : public GMVFile
{ // GMVFileAscii

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /** Constructor with name of GMV file.
   *
   * @param filename Name of GMV file
   */
  GMVFileAscii(const char* name);

  /// Default destructor 
  ~GMVFileAscii(void);

  /** Get header.
   *
   * @returns Header that appears in ASCII GMV file
   */
  static const char* header(void);

  /** Read ASCII GMV file.
   *
   * @coordinates Coordinates of vertices.
   * @param cells Indices of cell vertices.
   * @param materialIds Material identifiers for each cell.
   * @param meshDim Dimension of cells in mesh.
   * @param numVertices Number of vertices in mesh.
   * @param numCells Number of cells in mesh.
   * @param numCorners Number of vertices in each cell.
   */
  void read(scalar_array* coordinates,
	    int_array* cells,
	    int_array* materialIds,
	    int* meshDim,
	    int* spaceDim,
	    int* numVertices,
	    int* numCells,
	    int* numCorners);

  /** Write ASCII GMV file.
   *
   * @coordinates Coordinates of vertices.
   * @param cells Indices of cell vertices.
   * @param materialIds Material identifiers for each cell.
   * @param meshDim Dimension of cells in mesh.
   * @param spaceDim Number of coordinates per vertex.
   * @param numVertices Number of vertices in mesh.
   * @param numCells Number of cells in mesh.
   * @param numCorners Number of vertices in each cell.
   */
  void write(const scalar_array& coordinates,
	     const int_array& cells,
	     const int_array& materialIds,
	     const int meshDim,
	     const int spaceDim,
	     const int numVertices,
	     const int numCells,
	     const int numCorners);

// PRIVATE METHODS //////////////////////////////////////////////////////
 private :
  
  /** Read header.
   *
   * @param fin Input file stream
   */
  void _readHeader(std::ifstream& fin);

  /** Read vertices.
   *
   * @param fin Input file stream.
   * @param coordinates Coordinates of vertices.
   * @param numVertices Number of vertices.
   */
  void _readVertices(std::ifstream& fin,
		     scalar_array* coordinates,
		     int* numVertices,
		     int* spaceDim);

  /** Read cells.
   *
   * @param fin Input file stream
   * @param cells Indices of cell vertices.
   * @param numCells Number of cells in mesh.
   * @param numCorners Number of vertices in each cell.
   */
  void _readCells(std::ifstream& fin,
		  int_array* cells,
		  int* numCells,
		  int* numCorners);

  /** Read and discard variables associated with vertices.
   *
   * @param fin Input file stream
   * @param numVertices Number of vertices in mesh.
   * @param numCells Number of cells in mesh.
   */
  void _readVariables(std::ifstream& fin,
		      const int numVertices,
		      const int numCells);

  /** Read and discard material flags for vertices.
   *
   * @param fin Input file stream
   * @param numVertices Number of vertices in mesh.
   * @param numCells Number of cells in mesh.
   */
  void _readFlags(std::ifstream& fin,
		  const int numVertices,
		  const int numCells);

  /** Read material values for cells.
   *
   * @param fin Input file stream
   * @param materialIds Material identifiers for each cell.
   * @param numVertices Number of vertices in mesh.
   * @param numCells Number of cells in mesh.
   */
  void _readMaterials(std::ifstream& fin,
		      int_array* materialIds,
		      const int numVertices,
		      const int numCells);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :
  
  /** Header in ascii GMV file */
  static const char* _HEADER;

}; // GMVFileInAscii

#include "GMVFileAscii.icc" // inline methods

#endif // pylith_meshio_gmvfileascii


// End of file 
