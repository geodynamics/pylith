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

#if !defined(pylith_meshio_meshio_hh)
#define pylith_meshio_meshio_hh

#include <Mesh.hh> // PETSc Mesh

namespace pylith {
  namespace meshio {
    class MeshIO;
  } // meshio
} // pylith

class pylith::meshio::MeshIO
{ // MeshIO
// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :
  typedef ALE::Mesh        Mesh;
  typedef Mesh::sieve_type sieve_type;
  typedef Mesh::label_type label_type;
  
// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :
  typedef enum {VERTEX, CELL} PointType;

  /// Constructor
  MeshIO(void);

  /// Destructor
  virtual ~MeshIO(void);

  /** Set flag associated with building intermediate mesh topology
   *  elements.
   *
   * @param flag True to build intermediate topology, false not to build
   */
  void interpolate(const bool flag);

  /** Get flag associated with building intermediate mesh topology
   * elements.
   *
   * @returns True if building intermediate topology, false if not building
   */
  bool interpolate(void) const;

  /** Read mesh from file.
   *
   * @param mesh Pointer to PETSc mesh object
   */
  void read(ALE::Obj<Mesh>* mesh);

  /** Write mesh to file.
   *
   * @param mesh Pointer to PETSc mesh object
   */
  void write(ALE::Obj<Mesh>* mesh);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  /// Write mesh
  virtual
  void _write(void) const = 0;

  /// Read mesh
  virtual
  void _read(void) = 0;

  /** Get flag indicating whether indices start at 0 (True) or 1 (False).
   *
   * @returns True if indices start at 0, false if 1.
   */
  bool useIndexZero(void) const;

  /** Set flag indicating whether indices start at 0 (True) or 1 (False).
   *
   * @param flag True if indices start at 0, false if 1.
   */
  void useIndexZero(const bool flag);

  /** Get spatial dimension of mesh.
   *
   * @returns Spatial dimension of mesh
   */
  int getMeshDim(void) const;

  /** Build mesh topology and set vertex coordinates.
   *
   * @param coordinates Array of coordinates of vertices
   * @param numVertices Number of vertices
   * @param spaceDim Dimension of vector space for vertex coordinates
   * @param cells Array of indices of vertices in cells (first index is 0 or 1)
   * @param numCells Number of cells
   * @param numCorners Number of vertices per cell
   * @param meshDim Dimension of cells in mesh
   */
  void _buildMesh(const double* coordinates,
		  const int numVertices,
		  const int spaceDim,
		  const int* cells,
		  const int numCells,
		  const int numCorners,
		  const int meshDim);

  /** Get information about vertices in mesh.
   *
   * Method caller is responsible for memory management.
   *
   * @param pCoordinates Pointer to array of vertex coordinates
   * @param pNumVertices Pointer to number of vertices
   * @param pSpaceDim Poiner to dimension of vector space for coordinates
   */
  void _getVertices(double** pCoordinates,
		    int* pNumVertices,
		    int* pSpaceDim) const;

  /** Get information about cells in mesh.
   *
   * Method caller is responsible for memory management.
   *
   * @param pCells Pointer to array of indicates of vertices in each cell
   * @param pNumCells Pointer to number of cells in mesh
   * @param pNumCorners Pointer to number of vertices in each cell
   * @param pMeshDim Pointer to number of dimensions associated with cell
   */
  void _getCells(int** pCells,
		 int* pNumCells,
		 int* pNumCorners,
		 int* pMeshDim) const;

  /** Tag cells in mesh with material identifiers.
   *
   * @param materialIds Material identifiers [numCells]
   * @param numCells Number of cells
   */
  void _setMaterials(const int* materialIds,
		     const int numCells);

  /** Get material identifiers for cells.
   *
   * @param materialIds Material identifiers [numCells]
   * @param numCells Number of cells
   */
  void _getMaterials(int** pMaterialIds,
		     int* pNumCells) const;

  /** Build a point group
   *
   * @param name The group name
   * @param type The point type, e.g. VERTEX, CELL
   * @param numPoints The number of points
   * @param points An array of the points
   */
  void _buildGroup(const std::string& name,
                   const PointType type,
                   const int numPoints,
				   const int* points);

  /** Return all group names
   *
   */
  ALE::Obj<std::set<std::string> > _getGroups() const;

  /** Return a point group
   *
   * @param name The group name
   * @param type The point type, e.g. VERTEX, CELL
   * @param numPoints The number of points
   * @param points An array of the points
   */
  void _getGroup(const char *name,
                 PointType& type,
                 int& numPoints,
                 int *points[]) const;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  bool _useIndexZero; ///< Flag indicating if indicates start at 0 (T) or 1 (F)
  bool _interpolate; ///< True if building intermediate topology elements

  ALE::Obj<Mesh>* _mesh; ///< Pointer to PETSc mesh object

}; // MeshIO

#include "MeshIO.icc" // inline methods

#endif // pylith_meshio_meshio_hh

// End of file 
