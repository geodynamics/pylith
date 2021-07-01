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

#if !defined(pylith_meshio_meshdata1d_hh)
#define pylith_meshio_meshdata1d_hh

#include "MeshData.hh"

namespace pylith {
  namespace meshio {
     class MeshData1D;
  } // pylith
} // meshio

class pylith::meshio::MeshData1D : public MeshData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  MeshData1D(void);

  /// Destructor
  ~MeshData1D(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const int _numVertices; ///< Number of vertices
  static const int _spaceDim; ///< Number of dimensions in vertex coordinates
  static const int _numCells; ///< Number of cells
  static const int _cellDim; ///< Number of dimensions associated with cell
  static const int _numCorners; ///< Number of vertices in cell

  static const PylithScalar _vertices[]; ///< Pointer to coordinates of vertices
  static const int _cells[]; ///< Pointer to indices of vertices in cells
  static const int _materialIds[]; ///< Pointer to cell material identifiers

  static const int _groups[]; ///< Groups of points
  static const int _groupSizes[]; ///< Sizes of groups
  static const char* _groupNames[]; ///< Array of group names
  static const char* _groupTypes[]; ///< Array of group types
  static const int _numGroups; ///< Number of groups

  static const bool _useIndexZero; ///< First vertex is 0 if true, 1 if false

};

#endif // pylith_meshio_meshdata1d_hh

// End of file
