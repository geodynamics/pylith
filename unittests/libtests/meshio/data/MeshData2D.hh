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

#if !defined(pylith_meshio_meshdata2d_hh)
#define pylith_meshio_meshdata2d_hh

#include "MeshData.hh"

namespace pylith {
  namespace meshio {
     class MeshData2D;
  } // pylith
} // meshio

class pylith::meshio::MeshData2D : public MeshData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  MeshData2D(void);

  /// Destructor
  ~MeshData2D(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const int _numVertices; ///< Number of vertices
  static const int _spaceDim; ///< Number of dimensions in vertex coordinates
  static const int _numCells; ///< Number of cells
  static const int _cellDim; ///< Number of dimensions associated with cell
  static const int _numCorners; ///< Number of vertices in cell

  static const double _vertices[]; ///< Pointer to coordinates of vertices
  static const int _cells[]; ///< Pointer to indices of vertices in cells
  static const int _materialIds[]; ///< Pointer to cell material identifiers

  static const bool _useIndexZero; ///< First vertex is 0 if true, 1 if false

};

#endif // pylith_meshio_meshdata2d_hh

// End of file
