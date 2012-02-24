// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/MeshIOAscii.hh
 *
 * @brief C++ input/output manager for PyLith ASCII mesh files.
 */

#if !defined(pylith_meshio_meshioascii_hh)
#define pylith_meshio_meshioascii_hh

// Include directives ---------------------------------------------------
#include "MeshIO.hh" // ISA MeshIO

#include "spatialdata/utils/utilsfwd.hh" // USES LineParser

#include <iosfwd> // USES std::istream, std::ostream
#include <string> // HASA std::string

// MeshIOAscii ----------------------------------------------------------
/// C++ input/output manager for PyLith ASCII mesh files.
class pylith::meshio::MeshIOAscii : public MeshIO
{ // MeshIOAscii
  friend class TestMeshIOAscii; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  MeshIOAscii(void);

  /// Destructor
  ~MeshIOAscii(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
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

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /// Write mesh
  void _write(void) const;

  /// Read mesh
  void _read(void);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Read mesh vertices.
   *
   * @param parser Input parser.
   * @param coordinates Pointer to array of vertex coordinates
   * @param numVertices Pointer to number of vertices
   * @param spaceDim Pointer to dimension of coordinates vector space
   */
  void _readVertices(spatialdata::utils::LineParser& parser,
		     scalar_array* coordinates,
		     int* numVertices,
		     int* spaceDim) const;
  
  /** Write mesh vertices.
   *
   * @param fileout Output stream
   */
  void _writeVertices(std::ostream& fileout) const;
  
  /** Read mesh cells.
   *
   * @param parser Input parser.
   * @param pCells Pointer to array of indices of cell vertices
   * @param pMaterialIds Pointer to array of material identifiers
   * @param pNumCells Pointer to number of cells
   * @param pNumCorners Pointer to number of corners
   */
  void _readCells(spatialdata::utils::LineParser& parser,
		  int_array* pCells,
		  int_array* pMaterialIds,
		  int* numCells,
		  int* numCorners) const;
  
  /** Write mesh cells.
   *
   * @param fileout Output stream
   * @param cells Array of indices of cell vertices
   * @param numCells Number of cells
   * @param numCorners Number of corners
   */
  void _writeCells(std::ostream& fileout) const;
  
  /** Read a point group.
   *
   * @param parser Input parser.
   * @param mesh The mesh
   */
  void _readGroup(spatialdata::utils::LineParser& parser,
		  int_array* points,
                  GroupPtType* type,
                  std::string* name) const;
  
  /** Write a point group.
   *
   * @param fileout Output stream
   * @param name The group name
   */
  void _writeGroup(std::ostream& fileout,
		   const char* name) const;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::string _filename; ///< Name of file
  bool _useIndexZero; ///< Flag indicating if indicates start at 0 (T) or 1 (F)

  static
  const char *groupTypeNames[]; ///< Types of mesh groups.

}; // MeshIOAscii

#include "MeshIOAscii.icc" // inline methods

#endif // pylith_meshio_meshioascii_hh

// End of file 
