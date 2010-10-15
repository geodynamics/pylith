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

#include <portinfo>

#include "CellRefinerTri3.hh" // implementation of class methods

#include "MeshOrder.hh" // USES MeshOrder

#include <cassert> // USES assert()

#include <iostream> // TEMPORARY
// ----------------------------------------------------------------------
// Constructor
ALE::CellRefinerTri3::CellRefinerTri3(const mesh_type& mesh) :
 RefineEdges2(mesh)
{ // constructor
  assert(2 == mesh.getDimension());
} // constructor

// ----------------------------------------------------------------------
// Destructor
ALE::CellRefinerTri3::~CellRefinerTri3(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get number of refined cells for each original cell.
int
ALE::CellRefinerTri3::numNewCells(const point_type cell)
{ // numNewCells
  switch (_cellType(cell)) {
  case TRIANGLE:
    return 4;
  case LINE_COHESIVE_LAGRANGE:
    return 2;
  default:
    assert(0);
    throw ALE::Exception("Unknown cell type.");
  } // switch
} // numNewCells

// ----------------------------------------------------------------------
// Split cell into smaller cells of same type.
void
ALE::CellRefinerTri3::splitCell(const point_type cell,
				const point_type cone[],
				const int coneSize,
				point_type* curNewVertex)
{ // splitCell
  assert(curNewVertex);

  int numEdges = 0;
  const EdgeType* edges;
  
  switch (_cellType(cell)) {
  case TRIANGLE:
    _edges_TRIANGLE(&edges, &numEdges, cone, coneSize);
    break;
  case LINE_COHESIVE_LAGRANGE:
    _edges_LINE_COHESIVE_LAGRANGE(&edges, &numEdges, cone, coneSize);
    break;
  default:
    throw ALE::Exception("Unknown cell type.");
  } // switch

  for(int iEdge=0; iEdge < numEdges; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      // if vertex does not exist
      std::cout << "Edge: " << edges[iEdge] << ", new vertex: " << *curNewVertex << std::endl;
      _edgeToVertex[edges[iEdge]] = *curNewVertex;
      ++(*curNewVertex);
    } // if
  } // for
} // splitCell

// ----------------------------------------------------------------------
// Split cell into smaller cells of same type.
void
ALE::CellRefinerTri3::splitCellUncensored(const point_type cell,
					  const point_type cone[],
					  const int coneSize,
					  point_type* curNewVertex)
{ // splitCellUncensored
  assert(curNewVertex);

  int numEdges = 0;
  const EdgeType* edges;
  
  const bool uncensored = true;

  switch (_cellType(cell)) {
  case TRIANGLE:
    // No censored vertices on normal cell.
    break;
  case LINE_COHESIVE_LAGRANGE:
    _edges_LINE_COHESIVE_LAGRANGE(&edges, &numEdges, cone, coneSize, uncensored);
    break;
  default:
    throw ALE::Exception("Unknown cell type.");
  } // switch

  for(int iEdge=0; iEdge < numEdges; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      // if vertex does not exist
      std::cout << "Edge: " << edges[iEdge] << ", new vertex: " << *curNewVertex << std::endl;
      _edgeToVertex[edges[iEdge]] = *curNewVertex;
      ++(*curNewVertex);
    } // if
  } // for
} // splitCellUncensored

// ----------------------------------------------------------------------
// Get refined cells.
void
ALE::CellRefinerTri3::getNewCells(const point_type** cells,
				  int* numCells,
				  const point_type cell,
				  const point_type cone[],
				  const int coneSize,
				  const MeshOrder& orderOldMesh,
				  const MeshOrder& orderNewMesh)
{ // getNewCells
  assert(cells);
  assert(numCells);

  switch (_cellType(cell)) {
  case TRIANGLE: {
    const int coneVertexOffset = orderNewMesh.verticesNormal().min() - orderOldMesh.verticesNormal().min();
    _newCells_TRIANGLE(cells, numCells, cone, coneSize, coneVertexOffset);
    break;
  } // TRIANGLE
  case LINE_COHESIVE_LAGRANGE: {
    const int coneVertexOffsetNormal = orderNewMesh.verticesNormal().min() - orderOldMesh.verticesNormal().min();
    const int coneVertexOffsetCensored = orderNewMesh.verticesCensored().min() - orderOldMesh.verticesCensored().min();
    _newCells_LINE_COHESIVE_LAGRANGE(cells, numCells, cone, coneSize, coneVertexOffsetNormal, coneVertexOffsetCensored);
    break;
  } // LINE_COHESIVE_LAGRANGE
  default:
    throw ALE::Exception("Unknown cell type.");
  } // switch
} // getNewCells

// ----------------------------------------------------------------------
// Get cell type.
ALE::CellRefinerTri3::CellEnum
ALE::CellRefinerTri3::_cellType(const point_type cell)
{ // _cellType
  assert(!_mesh.getSieve().isNull());

  switch (_mesh.getSieve()->getConeSize(cell)) {
  case 3:
    return TRIANGLE;
  case 6:
    return LINE_COHESIVE_LAGRANGE;
  case 0: {
    std::ostringstream msg;
    std::cerr << "Internal error. Cone size for mesh point " << cell << " is zero. May be a vertex.";
    assert(0);
    throw ALE::Exception("Could not determine cell type during uniform global refinement.");
  } // case 0
  default : {
    std::ostringstream msg;
    std::cerr << "Internal error. Unknown cone size for mesh point " << cell << ". Unknown cell type.";
    assert(0);
    throw ALE::Exception("Could not determine cell type during uniform global refinement.");
  } // default
  } // switch
} // _cellType
  
// ----------------------------------------------------------------------
// Get edges of triangular cell.
void
ALE::CellRefinerTri3::_edges_TRIANGLE(const EdgeType** edges,
				      int* numEdges,
				      const point_type cone[],
				      const int coneSize)
{ // _edges_TRIANGLE
  static EdgeType triEdges[3];
  
  assert(coneSize == 3);
  triEdges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
  triEdges[1] = EdgeType(std::min(cone[1], cone[2]), std::max(cone[1], cone[2]));
  triEdges[2] = EdgeType(std::min(cone[2], cone[0]), std::max(cone[2], cone[0]));
  *numEdges = 3;
  *edges    = triEdges;
} // _edges_TRIANGLE
  
// ----------------------------------------------------------------------
// Get edges of line cohesive cell with Lagrange multipler vertices.
void
ALE::CellRefinerTri3::_edges_LINE_COHESIVE_LAGRANGE(const EdgeType** edges,
						    int* numEdges,
						    const point_type cone[],
						    const int coneSize,
						    const bool uncensored)
{ // _edges_LINE_COHESIVE_LAGRANGE
  if (uncensored) {
    // Include all edges
    static EdgeType lineEdges[3];

    assert(coneSize == 6);
    lineEdges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
    lineEdges[1] = EdgeType(std::min(cone[2], cone[3]), std::max(cone[2], cone[3]));
    lineEdges[2] = EdgeType(std::min(cone[4], cone[5]), std::max(cone[4], cone[5]));
    *numEdges = 3;
    *edges    = lineEdges;
  } else {
    // Omit edges with censored (Lagrange multiplier) vertices.
    static EdgeType lineEdges[2];

    assert(coneSize == 6);
    lineEdges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
    lineEdges[1] = EdgeType(std::min(cone[2], cone[3]), std::max(cone[2], cone[3]));
    *numEdges = 2;
    *edges    = lineEdges;
  } // if/else
} // _edges_LINE_COHESIVE_LAGRANGE
  
// ----------------------------------------------------------------------
// Get new cells from refinement of a triangular cell.
void
ALE::CellRefinerTri3::_newCells_TRIANGLE(const point_type** cells,
					 int *numCells,
					 const point_type cone[],
					 const int coneSize,
					 const int coneVertexOffset)
{ // _newCells_TRIANGLE
  const int coneSizeTri3 = 3;
  const int numEdgesTri3 = 3;
  const int numNewCells = 4;
  const int numNewVertices = 3;

  int numEdges = 0;
  const EdgeType  *edges;
  _edges_TRIANGLE(&edges, &numEdges, cone, coneSize);
  assert(numEdgesTri3 == numEdges);

  static point_type triCells[numNewCells*coneSizeTri3];
  point_type newVertices[numNewVertices];
  for(int iEdge=0, iNewVertex=0; iEdge < numEdgesTri3; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      throw ALE::Exception("Missing edge in refined mesh");
    } // if
    newVertices[iNewVertex++] = _edgeToVertex[edges[iEdge]];
  } // for

  // new cell 0
  triCells[0*3+0] = cone[0] + coneVertexOffset;
  triCells[0*3+1] = newVertices[0];
  triCells[0*3+2] = newVertices[2];

  // new cell 1
  triCells[1*3+0] = newVertices[0];
  triCells[1*3+1] = newVertices[1];
  triCells[1*3+2] = newVertices[2];

  // new cell 2
  triCells[2*3+0] = cone[1] + coneVertexOffset;
  triCells[2*3+1] = newVertices[1];
  triCells[2*3+2] = newVertices[0];

  // new cell 3
  triCells[3*3+0] = cone[2] + coneVertexOffset;
  triCells[3*3+1] = newVertices[2];
  triCells[3*3+2] = newVertices[1];

  *numCells = numNewCells;
  *cells    = triCells;
} // _newCells_TRIANGLE
  
// ----------------------------------------------------------------------
// Get new cells from refinement of a line cohseive cell with Lagrange
// multiplier vertices.
void
ALE::CellRefinerTri3::_newCells_LINE_COHESIVE_LAGRANGE(const point_type** cells,
						       int *numCells,
						       const point_type cone[],
						       const int coneSize,
						       const int coneVertexOffsetNormal,
						       const int coneVertexOffsetCensored)
{ // _newCells_LINE_COHESIVE_LAGRANGE
  const int coneSizeLine6 = 6;
  const int numEdgesLine6 = 3;
  const int numNewCells = 2;
  const int numNewVertices = 3;

  int numEdges = 0;
  const EdgeType *edges;
  _edges_LINE_COHESIVE_LAGRANGE(&edges, &numEdges, cone, coneSize, true);
  assert(numEdgesLine6 == numEdges);

  static point_type lineCells[numNewCells*coneSizeLine6];
  point_type newVertices[numNewVertices];
  for(int iEdge=0, iNewVertex=0; iEdge < numEdgesLine6; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      throw ALE::Exception("Missing edge in refined mesh");
    } // if
    newVertices[iNewVertex++] = _edgeToVertex[edges[iEdge]];
  } // for

  // new cell 0
  lineCells[0*6+0] = cone[0] + coneVertexOffsetNormal;
  lineCells[0*6+1] = newVertices[0];
  lineCells[0*6+2] = cone[2] + coneVertexOffsetNormal;
  lineCells[0*6+3] = newVertices[1];
  lineCells[0*6+4] = cone[4] + coneVertexOffsetCensored;
  lineCells[0*6+5] = newVertices[2];

  // new cell 1
  lineCells[1*6+0] = newVertices[0];
  lineCells[1*6+1] = cone[1] + coneVertexOffsetNormal;
  lineCells[1*6+2] = newVertices[1];
  lineCells[1*6+3] = cone[3] + coneVertexOffsetNormal;
  lineCells[1*6+4] = newVertices[2];
  lineCells[1*6+5] = cone[5] + coneVertexOffsetCensored;
  
  *numCells = 2;
  *cells    = lineCells;
} // _newCells_LINE_COHESIVE_LAGRANGE


// End of file 
