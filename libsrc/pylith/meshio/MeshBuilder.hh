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

/**
 * @file libsrc/meshio/MeshBuilder.hh
 *
 * @brief Helper class for constructing PETSc mesh data structures.
 */

#if !defined(pylith_meshio_meshbuilder_hh)
#define pylith_meshio_meshbuilder_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "pylith/utils/arrayfwd.hh" // USES scalar_array, int_array,
                                    // string_vector
#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

#include "pylith/topology/Mesh.hh" // USES Mesh

// MeshBuilder ----------------------------------------------------------
/// Helper class for constructing PETSc mesh data structures.
class pylith::meshio::MeshBuilder
{ // MeshBuilder

// PUBLIC TYPEDEFS ////////////////////////////////////////////////////
public:

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

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
 * @param isParallel Create parallel mesh if true, otherwise only build
 *   mesh on proc 0.
 */
static
void buildMesh(topology::Mesh* mesh,
               scalar_array* coordinates,
               const int numVertices,
               int spaceDim,
               const int_array& cells,
               const int numCells,
               const int numCorners,
               const int meshDim,
               const bool isParallel =false);
}; // MeshBuilder

#endif // pylith_meshio_meshbuilder_hh

// End of file
