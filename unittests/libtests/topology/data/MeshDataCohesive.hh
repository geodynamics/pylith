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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_topology_meshdatacohesive_hh)
#define pylith_topology_meshdatacohesive_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
  namespace topology {
     class MeshDataCohesive;
  } // pylith
} // topology

class pylith::topology::MeshDataCohesive
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  MeshDataCohesive(void);

  /// Destructor
  ~MeshDataCohesive(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  // Input information
  char* filename; ///< Filename of mesh file.
  int refineLevel; ///< Refinement level.
  char* faultA; ///< Vertex group associated with fault A (0 if no fault).
  char* faultB; ///< Vertex group associated with fault B (0 if no fault).

  // Output information
  int numVertices; ///< Number of vertices
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int cellDim; ///< Number of dimensions associated with cell.
  int numCells; ///< Number of cells
  int numCorners; ///< Number of vertices in cell.
  int numCellsCohesive; ///< Number of cohesive cells.
  int numCornersCohesive; ///< Number of vertices in cohesive cell.

  PylithScalar* vertices; ///< Pointer to coordinates of vertices
  int* cells; ///< Pointer to indices of vertices in cells
  int* cellsCohesive; ///< Pointer to indices of vertices in cells
  int* materialIds; ///< Pointer to cell material identifiers

  int* groups; ///< Array of pointers to indices of points in groups
  int* groupSizes; ///< Array of sizes of each group
  char** groupNames; ///< Array of group names
  char** groupTypes; ///< Array of group types
  int numGroups; ///< Number of groups

};

#endif // pylith_topology_meshdatacohesive_hh

// End of file
