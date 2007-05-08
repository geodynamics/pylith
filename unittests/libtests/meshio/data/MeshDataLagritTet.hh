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

#if !defined(pylith_meshio_meshdatalagrittet_hh)
#define pylith_meshio_meshdatalagrittet_hh

#include "MeshData.hh"

namespace pylith {
  namespace meshio {
     class MeshDataLagritTet;
  } // pylith
} // meshio

class pylith::meshio::MeshDataLagritTet : public MeshData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  MeshDataLagritTet(void);

  /// Destructor
  ~MeshDataLagritTet(void);

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

  static const int _groups[]; ///< Groups of points
  static const int _groupSizes[]; ///< Sizes of groups
  static const char* _groupNames[]; ///< Array of group names
  static const char* _groupTypes[]; ///< Array of group types
  static const int _numGroups; ///< Number of groups

  static const bool _useIndexZero; ///< First vertex is 0 if true, 1 if false

};

#endif // pylith_meshio_meshdatalagrittet_hh

// End of file
