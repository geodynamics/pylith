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

/**
 * @file pylith/meshio/MeshIOCubit.hh
 *
 * @brief C++ input/output manager for CUBIT Exodus II files.
 */

#if !defined(pylith_meshio_meshiocubit_hh)
#define pylith_meshio_meshiocubit_hh

// Include directives ---------------------------------------------------
#include "MeshIO.hh" // ISA MeshIO

#include <string> // HASA std::string

// Forward declarations -------------------------------------------------
class NcFile; // netcdf file

// MeshIOCubit ----------------------------------------------------------
class pylith::meshio::MeshIOCubit : public MeshIO
{ // MeshIOCubit
  friend class TestMeshIOCubit; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  MeshIOCubit(void);

  /// Destructor
  ~MeshIOCubit(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Set filename for Cubit file.
   *
   * @param filename Name of file
   */
  void filename(const char* name);

  /** Get filename of Cubit file.
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
   * @param ncfile Cubit Exodus file.
   * @param coordinates Pointer to array of vertex coordinates.
   * @param numVertices Pointer to number of vertices.
   * @param spaceDim Pointer to dimension of coordinates vector space.
   */
  void _readVertices(NcFile& filein,
		     double_array* coordinates,
		     int* numVertices,
		     int* spaceDim) const;
  
  /** Read mesh cells.
   *
   * @param ncfile Cubit Exodus file.
   * @param pCells Pointer to array of indices of cell vertices
   * @param pMaterialIds Pointer to array of material identifiers
   * @param pNumCells Pointer to number of cells
   * @param pNumCorners Pointer to number of corners
   */
  void _readCells(NcFile& filein,
		  int_array* pCells,
		  int_array* pMaterialIds,
		  int* numCells,
		  int* numCorners) const;
  
  /** Read point groups.
   *
   * @param ncfile Cubit Exodus file.
   */
  void _readGroups(NcFile& filein);
  
  /** Write mesh dimensions.
   *
   * @param ncfile Cubit Exodus file.
   */
  void _writeDimensions(NcFile& ncfile) const;
  
  /** Write mesh variables.
   *
   * @param ncfile Cubit Exodus file.
   */
  void _writeVariables(NcFile& ncfile) const;
  
  /** Write mesh attributes.
   *
   * @param ncfile Cubit Exodus file.
   */
  void _writeAttributes(NcFile& ncfile) const;

  /** Reorder vertices in cells to match PyLith conventions.
   *
   * @param cells Array of vertex indices for each cell [numCells*numCorners].
   * @param numCells Number of cells.
   * @param numCorners Number of vertices per cell.
   * @param meshDim Spatial dimension of mesh.
   */
  static
  void _orientCells(int_array* const cells,
		    const int numCells,
		    const int numCorners,
		    const int meshDim);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::string _filename; ///< Name of file

}; // MeshIOCubit

#include "MeshIOCubit.icc" // inline methods

#endif // pylith_meshio_meshiocubit_hh

// End of file 
