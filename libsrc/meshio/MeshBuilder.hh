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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/MeshBuilder.hh
 *
 * @brief Helper class for constructing Sieve mesh data structures.
 */

#if !defined(pylith_meshio_meshbuilder_hh)
#define pylith_meshio_meshbuilder_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "pylith/utils/arrayfwd.hh" // USES double_array, int_array,
                                    // string_vector
#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/sievetypes.hh" // USE SieveFlexMesh

// MeshBuilder ----------------------------------------------------------
/// Helper class for constructing Sieve mesh data structures.
class pylith::meshio::MeshBuilder
{ // MeshBuilder

  // PUBLIC TYPEDEFS ////////////////////////////////////////////////////
public :

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Build mesh topology and set vertex coordinates.
   *
   * All mesh information must use zero based indices. In other words,
   * the lowest index MUST be 0 not 1.
   *
   * @param mesh PyLith finite-element mesh.
   * @param coordinates Array of coordinates of vertices.
   * @param numVertices Number of vertices.
   * @param spaceDim Dimension of vector space for vertex coordinates.
   * @param cells Array of indices of vertices in cells (first index is 0).
   * @param numCells Number of cells.
   * @param numCorners Number of vertices per cell.
   * @param meshDim Dimension of cells in mesh.
   */
  static
  void buildMesh(topology::Mesh* mesh,
		 double_array* coordinates,
		 const int numVertices,
		 const int spaceDim,
		 const int_array& cells,
		 const int numCells,
		 const int numCorners,
		 const int meshDim,
		 const bool interpolate);

  /** Build fault mesh topology.
   *
   * All mesh information must use zero based indices. In other words,
   * the lowest index MUST be 0 not 1.
   *
   * @param fault PETSc mesh for fault.
   * @param faultBd PETSc mesh for fault boundary.
   * @param coordinates Array of coordinates of vertices.
   * @param numVertices Number of vertices.
   * @param spaceDim Dimension of vector space for vertex coordinates.
   * @param cells Array of indices of vertices in cells (first index is 0).
   * @param numCells Number of cells.
   * @param numCorners Number of vertices per cell.
   * @param firstCell Label of first cell.
   * @param cells Array of indices of vertices for fault surface cells
   * (first index is 0).
   * @param meshDim Dimension of cells in mesh.
   */
  static
  void buildFaultMesh(const ALE::Obj<topology::Mesh::SieveMesh>& fault,
		      ALE::Obj<SieveFlexMesh>& faultBd,
		      const double_array& coordinates,
		      const int numVertices,
		      const int spaceDim,
		      const int_array& cells,
		      const int numCells,
		      const int numCorners,
		      const int firstCell,
		      const int_array& faceCells,
		      const int meshDim);


}; // MeshBuilder

#endif // pylith_meshio_meshbuilder_hh

// End of file 


