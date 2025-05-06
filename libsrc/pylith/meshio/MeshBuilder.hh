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
#include "pylith/utils/arrayfwd.hh" // HASA scalar_array, int_array,  string_vector

#include "pylith/topology/Mesh.hh" // USES Mesh

/// Helper class for constructing PETSc mesh data structures.
class pylith::meshio::MeshBuilder {
    // PUBLIC ENUMS /////////////////////////////////////////////////////
public:

    enum shape_t {
        POINT=0,
        LINE=1,
        TRIANGLE=2,
        QUADRILATERAL=3,
        TETRAHEDRON=4,
        HEXAHEDRON=5,
        SHELL4=6, // Exodus-II
    };

    struct Topology {
        size_t dimension;
        size_t numCells;
        size_t numCorners;
        shape_t cellShape;
        int_array cells;

        Topology(void) :
            dimension(0),
            numCells(0),
            numCorners(0),
            cellShape(POINT) {}


        Topology(const size_t _dimension,
                 const size_t _numCells,
                 const size_t _numCorners,
                 const shape_t _cellShape,
                 const int* _cells) :
            dimension(_dimension),
            numCells(_numCells),
            numCorners(_numCorners),
            cellShape(_cellShape),
            cells(_cells, _numCells*_numCorners) {}


    };

    struct Geometry {
        size_t numVertices;
        size_t spaceDim;
        scalar_array vertices;

        Geometry(void) :
            numVertices(0),
            spaceDim(0) {}


        Geometry(const size_t _numVertices,
                 const size_t _spaceDim,
                 const PylithScalar* _vertices) :
            numVertices(_numVertices),
            spaceDim(_spaceDim),
            vertices(_vertices, _numVertices*_spaceDim) {}


    };

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Build mesh topology and set vertex coordinates.
     *
     * All mesh information must use zero based indices. In other words,
     * the lowest index MUST be 0 not 1.
     *
     * @param[inout] mesh PyLith finite-element mesh.
     * @param[in] topology Mesh topology information.
     * @param[in] geometry Mesh geometry information.
     * @param[in] isParallel Create parallel mesh if true, otherwise only build
     *   mesh on proc 0.
     */
    static
    void buildMesh(pylith::topology::Mesh* mesh,
                   const Topology& topology,
                   const Geometry& geometry,
                   const bool isParallel=false);

    /** Tag cells in mesh with material identifiers.
     *
     * @param[inout] mesh Finite-element mesh.
     * @param[in] materialIds Material identifiers [numCells]
     */
    static
    void setMaterials(pylith::topology::Mesh* mesh,
                      const int_array& materialIds);

    /** Build a point group of vertices.
     *
     * The indices in the points array must use zero based indices. In
     * other words, the lowest index MUST be 0 not 1.
     *
     * @param[inout] mesh PyLith finite-element mesh.
     * @param[in] name The group name
     * @param[in] points An array of the points in the group.
     * @param[in] labelValue Value of group label in mesh.
     */
    static
    void setVertexGroup(pylith::topology::Mesh* mesh,
                        const char* name,
                        const int_array& points,
                        const int labelValue=1);

    /** Build a point group of faces from cell+vertices.
     *
     * The indices in the points array must use zero based indices. In
     * other words, the lowest index MUST be 0 not 1.
     *
     * @param[inout] mesh PyLith finite-element mesh.
     * @param[in] name The group name
     * @param[in] points An array of the points in the group.
     * @param[in] faceShape Shape of face.
     * @param[in] labelValue Value of group label in mesh.
     */
    static
    void setFaceGroupFromCellVertices(pylith::topology::Mesh* mesh,
                                      const char* name,
                                      const int_array& points,
                                      const shape_t faceShape,
                                      const int labelValue=1);

    /** Build a point group of faces from cell+side.
     *
     * The indices in the points array must use zero based indices. In
     * other words, the lowest index MUST be 0 not 1.
     *
     * @param[inout] mesh PyLith finite-element mesh.
     * @param[in] name The group name
     * @param[in] points An array of the points in the group.
     * @param[in] labelValue Value of group label in mesh.
     */
    static
    void setFaceGroupFromCellSide(pylith::topology::Mesh* mesh,
                                  const char* name,
                                  const int_array& points,
                                  const int labelValue=1);

    /** Get information about vertices in mesh.
     *
     * @param[inout] geometry Mesh geometry information.
     * @param[in] mesh Finite-element mesh.
     */
    static
    void getVertices(Geometry* geometry,
                     const pylith::topology::Mesh& mesh);

    /** Get information about cells in mesh.
     *
     * The indices of the vertices will be returned using zero based
     * indices. In other words, the first vertex will be referred to
     * using an index of 0.
     *
     * @param[inout] meshData Mesh information.
     * @param[in] mesh Finite-element mesh.
     */
    static
    void getCells(Topology* topology,
                  const pylith::topology::Mesh& mesh);

    /** Get material identifiers for cells.
     *
     * @param[out] materialIds Material identifiers [numCells]
     * @param[in] mesh Finite-element mesh.
     */
    static
    void getMaterials(int_array* materialIds,
                      const pylith::topology::Mesh& mesh);

    /** Get names of vertex groups in mesh.
     *
     * @returns Array of group names.
     * @param[in] mesh Finite-element mesh.
     */
    static
    void getVertexGroupNames(string_vector* names,
                             const pylith::topology::Mesh& mesh);

    /** Get names of face groups in mesh.
     *
     * @returns Array of group names.
     * @param[in] mesh Finite-element mesh.
     */
    static
    void getFaceGroupNames(string_vector* names,
                           const pylith::topology::Mesh& mesh);

    /** Get a point group of vertices.
     *
     * The indices in the points array use zero based indices and have the
     * offset removed, so the first vertex is the mesh is 0.
     *
     * @param[in] points An array of the points in the group.
     * @param[inout] mesh PyLith finite-element mesh.
     * @param[in] name The group name
     * @param[in] labelValue Value of group label in mesh.
     */
    static
    void getVertexGroup(int_array* points,
                        const pylith::topology::Mesh& mesh,
                        const char* name,
                        const int labelValue=1);

    /** Get a point group of faces.
     *
     * The indices in the points array use zero based indices and have the
     * offset removed, so the first vertex is the mesh is 0.
     *
     * @param[in] points An array of the points in the group.
     * @param[inout] mesh PyLith finite-element mesh.
     * @param[in] name The group name
     * @param[in] labelValue Value of group label in mesh.
     */
    static
    void getFaceGroup(int_array* points,
                      const pylith::topology::Mesh& mesh,
                      const char* name,
                      const int labelValue=1);

    /* Get cell shape from dimension and number of corners.
     *
     * @returns Cell shape.
     * @param[in] cellDim Cell dimension.
     * @param[in] numCorners Number of corners (vertices) in cell.
     */
    static
    shape_t cellShapeFromCorners(const size_t cellDim,
                                 const size_t numCorners);

    /** Get face shape from cell shape.
     *
     * @returns Shape of cell face.
     * @param[in] cellShape Cell shape.
     */
    static
    shape_t faceShapeFromCellShape(const shape_t cellShape);

    /** Get number of face vertices given face shape.
     *
     * @returns Number of face vertices.
     * @param[in] faceShape Face shape.
     */
    static
    size_t getNumVerticesFace(const pylith::meshio::MeshBuilder::shape_t faceShape);

}; // MeshBuilder

// End of file
