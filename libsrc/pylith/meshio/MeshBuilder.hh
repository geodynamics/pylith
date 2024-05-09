// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/meshio/meshiofwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "pylith/utils/arrayfwd.hh" // USES scalar_array, int_array,
                                    // string_vector
#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

#include "pylith/topology/Mesh.hh" // USES Mesh

/// Helper class for constructing PETSc mesh data structures.
class pylith::meshio::MeshBuilder {
    // PUBLIC ENUMS /////////////////////////////////////////////////////
public:

    /// Type of points in a group.
    enum GroupPtType {
        VERTEX=0,
        FACE=1,
        CELL=2,
    }; // GroupPtType

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Build mesh topology and set vertex coordinates.
     *
     * All mesh information must use zero based indices. In other words,
     * the lowest index MUST be 0 not 1.
     *
     * @param[inout] mesh PyLith finite-element mesh.
     * @param[in] coordinates Array of coordinates of vertices.
     * @param[in] numVertices Number of vertices.
     * @param[in] spaceDim Dimension of vector space for vertex coordinates.
     * @param[in] cells Array of indices of vertices in cells (first index is 0).
     * @param[in] numCells Number of cells.
     * @param[in] numCorners Number of vertices per cell.
     * @param[in] meshDim Dimension of cells in mesh.
     * @param[in] isParallel Create parallel mesh if true, otherwise only build
     *   mesh on proc 0.
     */
    static
    void buildMesh(pylith::topology::Mesh* mesh,
                   scalar_array* coordinates,
                   const int numVertices,
                   int spaceDim,
                   const int_array& cells,
                   const int numCells,
                   const int numCorners,
                   const int meshDim,
                   const bool isParallel=false);

    /** Build a point group
     *
     * The indices in the points array must use zero based indices. In
     * other words, the lowest index MUST be 0 not 1.
     *
     * @param[inout] mesh PyLith finite-element mesh.
     * @param[in] name The group name
     * @param[in] groupType The point type, e.g. VERTEX, CELL
     * @param[in] points An array of the points in the group.
     */
    static
    void setGroup(pylith::topology::Mesh* mesh,
                  const char* name,
                  const GroupPtType groupType,
                  const int_array& points);

    /** Get a point group as an array of points.
     *
     * The indices in the points array must use zero based indices. In
     * other words, the lowest index MUST be 0 not 1.
     *
     * @param[inout] mesh PyLith finite-element mesh.
     * @param[in] name The group name
     * @param[in] groupType The point type, e.g. VERTEX, CELL
     * @param[in] points An array of the points in the group.
     */
    static
    void getGroup(GroupPtType* groupType,
                  int_array* points,
                  const pylith::topology::Mesh& mesh,
                  const char* name);

}; // MeshBuilder

// End of file
