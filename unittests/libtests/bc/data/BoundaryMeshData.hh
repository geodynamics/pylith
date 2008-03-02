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

#if !defined(pylith_bc_boundarymeshdata_hh)
#define pylith_bc_boundarymeshdata_hh

namespace pylith {
  namespace bc {
     class BoundaryMeshData;
  } // pylith
} // bc

class pylith::bc::BoundaryMeshData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  BoundaryMeshData(void);

  /// Destructor
  ~BoundaryMeshData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* filename; ///< Name of file with input mesh.
  char* bcLabel; ///< Name of group of vertices for bc.
  char* faultLabel; ///< Name of group of vertices for fault.
  int faultId; ///< Material identifier for fault.

  int numCorners; ///< Number of vertices in cells of boundary mesh.
  int numCells; ///< Number of cells in boundary mesh.

  /// @name Boundary mesh without fault.
  //@{
  int numVerticesNoFault; ///< Number of vertices.
  int* verticesNoFault; ///< Array of vertex labels.
  int* cellsNoFault; ///< Array of vertex labels for cells.
  //@}

  /// @name Boundary mesh without fault.
  //@{
  int numVerticesFault; ///< Number of vertices.
  int* verticesFault; ///< Array of vertex labels.
  int* cellsFault; ///< Array of vertex labels for cells.
  //@}

};

#endif // pylith_bc_boundarymeshdata_hh


// End of file
