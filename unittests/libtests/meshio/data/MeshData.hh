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

#if !defined(pylith_meshio_meshdata_hh)
#define pylith_meshio_meshdata_hh

namespace pylith {
  namespace meshio {
     class MeshData;
  } // pylith
} // meshio

class pylith::meshio::MeshData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  MeshData(void);

  /// Destructor
  ~MeshData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  int numVertices; ///< Number of vertices
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int numCells; ///< Number of cells
  int cellDim; ///< Number of dimensions associated with cell
  int numCorners; ///< Number of vertices in cell

  double* vertices; ///< Pointer to coordinates of vertices
  int* cells; ///< Pointer to indices of vertices in cells

  bool useIndexZero; ///< Indices start with 0 if true, 1 if false

  // :TODO:
  // Add groups of vertices

};

#endif // pylith_meshio_meshdata_hh

// End of file
