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

#include "pylith/utils/arrayfwd.hh" // USES double_array, int_array,
                                    // string_vector

#include "pylith/utils/sievetypes.hh" // USES Obj, PETSc Mesh

#include "Mesh.hh" // Needs ALE::Mesh

namespace pylith {
  namespace meshio {
    class MeshIO;
  } // meshio
} // pylith

class pylith::meshio::MeshIO
{ // MeshIO

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  /// Type of points in a group.
  typedef enum { VERTEX, CELL } GroupPtType;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :
  /// Constructor
  MeshIO(void);

  /// Destructor
  virtual ~MeshIO(void);

  /** Set debug flag for mesh.
   *
   * @param flag True to print debugging information.
   */
  void debug(const bool flag);

  /** Get debug flag for mesh.
   *
   * @returns True if debugging is on.
   */
  bool debug(void) const;

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
  void read(ALE::Obj<pylith::Mesh>* mesh);

  /** Write mesh to file.
   *
   * @param mesh Pointer to PETSc mesh object
   */
  void write(ALE::Obj<pylith::Mesh>* mesh);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  /// Write mesh
  virtual
  void _write(void) const = 0;

  /// Read mesh
  virtual
  void _read(void) = 0;

  /** Get spatial dimension of mesh.
   *
   * @returns Spatial dimension of mesh
   */
  int getMeshDim(void) const;

  /** Build mesh topology and set vertex coordinates.
   *
   * All mesh information must use zero based indices. In other words,
   * the lowest index MUST be 0 not 1.
   *
   * @param coordinates Array of coordinates of vertices
   * @param numVertices Number of vertices
   * @param spaceDim Dimension of vector space for vertex coordinates
   * @param cells Array of indices of vertices in cells (first index is 0 or 1)
   * @param numCells Number of cells
   * @param numCorners Number of vertices per cell
   * @param meshDim Dimension of cells in mesh
   */
  void _buildMesh(const double_array& coordinates,
		  const int numVertices,
		  const int spaceDim,
		  const int_array& cells,
		  const int numCells,
		  const int numCorners,
		  const int meshDim);

  static void
  _buildFaultMesh(const double_array& coordinates,
                  const int numVertices,
                  const int spaceDim,
                  const int_array& cells,
                  const int numCells,
                  const int numCorners,
                  const int firstCell,
                  const int_array& faceCells,
                  const int meshDim,
                  const Obj<Mesh>& fault,
                  Obj<ALE::Mesh>& faultBd);

  /** Get information about vertices in mesh.
   *
   * Method caller is responsible for memory management.
   *
   * @param coordinates Pointer to array of vertex coordinates
   * @param numVertices Pointer to number of vertices
   * @param spaceDim Poiner to dimension of vector space for coordinates
   */
  void _getVertices(double_array* coordinates,
		    int* numVertices,
		    int* spaceDim) const;

  /** Get information about cells in mesh.
   *
   * The indices of the vertices will be returned using zero based
   * indices. In other words, the first vertex will be referred to
   * using an index of 0.
   *
   * Method caller is responsible for memory management.
   *
   * @param cells Pointer to array of indicates of vertices in each cell
   * @param numCells Pointer to number of cells in mesh
   * @param numCorners Pointer to number of vertices in each cell
   * @param meshDim Pointer to number of dimensions associated with cell
   */
  void _getCells(int_array* cells,
		 int* numCells,
		 int* numCorners,
		 int* meshDim) const;

  /** Tag cells in mesh with material identifiers.
   *
   * @param materialIds Material identifiers [numCells]
   */
  void _setMaterials(const int_array& materialIds);

  /** Get material identifiers for cells.
   *
   * @param materialIds Material identifiers [numCells]
   */
  void _getMaterials(int_array* pMaterialIds) const;

  /** Build a point group
   *
   * The indices in the points array must use zero based indices. In
   * other words, the lowest index MUST be 0 not 1.
   *
   * @param name The group name
   * @param type The point type, e.g. VERTEX, CELL
   * @param points An array of the points in the group.
   */
  void _setGroup(const std::string& name,
		 const GroupPtType type,
		 const int_array& points);

  /** Get names of all groups in mesh.
   *
   * @returns Array of group names.
   */
  void _getGroupNames(string_vector* names) const;

  /** Return a point group
   *
   * The indices of the vertices will be returned using zero based
   * indices. In other words, the first vertex will be referred to
   * using an index of 0.
   *
   * @param points An array of the points in the group
   * @param type The point type, e.g. VERTEX, CELL
   * @param name The group name
   */
  void _getGroup(int_array* points,
		 GroupPtType* type,
		 const char *name) const;

  /** Create empty groups on other processes
   */
  void _distributeGroups();

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  bool _debug; ///< True to turn of mesh debugging output
  bool _interpolate; ///< True if building intermediate topology elements

  ALE::Obj<pylith::Mesh>* _mesh; ///< Pointer to PETSc mesh object

}; // MeshIO

#include "MeshIO.icc" // inline methods

#endif // pylith_meshio_meshio_hh

// End of file 
