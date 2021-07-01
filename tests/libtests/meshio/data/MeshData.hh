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

#if !defined(pylith_meshio_meshdata_hh)
#define pylith_meshio_meshdata_hh

#include "pylith/utils/types.hh" // USES PylithScalar

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

  PylithScalar* vertices; ///< Pointer to coordinates of vertices
  int* cells; ///< Pointer to indices of vertices in cells
  int* materialIds; ///< Pointer to cell material identifiers

  int* groups; ///< Array of pointers to indices of points in groups
  int* groupSizes; ///< Array of sizes of each group
  char** groupNames; ///< Array of group names
  char** groupTypes; ///< Array of group types
  int numGroups; ///< Number of groups

  bool useIndexZero; ///< Indices start with 0 if true, 1 if false

};

#endif // pylith_meshio_meshdata_hh

// End of file
