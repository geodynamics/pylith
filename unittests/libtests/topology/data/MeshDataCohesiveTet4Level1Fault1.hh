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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_topology_meshdatacohesivetet4level1fault1_hh)
#define pylith_topology_meshdatacohesivetet4level1fault1_hh

#include "MeshDataCohesive.hh"

namespace pylith {
  namespace topology {
     class MeshDataCohesiveTet4Level1Fault1;
  } // pylith
} // topology

class pylith::topology::MeshDataCohesiveTet4Level1Fault1 : public MeshDataCohesive
{ // MeshDataCohesiveTet4Level1Fault1

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  MeshDataCohesiveTet4Level1Fault1(void);

  /// Destructor
  ~MeshDataCohesiveTet4Level1Fault1(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const char* _filename; ///< Filename of mesh file.
  static const int _refineLevel; ///< Refinement level.
  static const char* _faultA; ///< Vertex group associated with fault A (0 if no fault).
  static const char* _faultB; ///< Vertex group associated with fault B (0 if no fault).

  static const int _numVertices; ///< Number of vertices
  static const int _spaceDim; ///< Number of dimensions in vertex coordinates
  static const int _cellDim; ///< Number of dimensions associated with cell
  static const int _numCells; ///< Number of cells
  static const int _numCorners; ///< Number of vertices in cell
  static const int _numCellsCohesive; ///< Number of cohesive cells.
  static const int _numCornersCohesive; ///< Number of vertices in cohesive cell.

  static const int _matIdSum; ///< Sum of material identifiers (checksum).

  static const int _groupSizes[]; ///< Sizes of groups
  static const char* _groupNames[]; ///< Array of group names
  static const char* _groupTypes[]; ///< Array of group types
  static const int _numGroups; ///< Number of groups

}; // MeshDataCohesiveTet4Level1Fault1

#endif // pylith_topology_meshdatacohesivetet4level1fault1_hh

// End of file
