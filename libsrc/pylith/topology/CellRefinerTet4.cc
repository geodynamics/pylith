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

#include <portinfo>

#include "CellRefinerTet4.hh" // implementation of class methods

#include "MeshOrder.hh" // USES MeshOrder

#include <cassert> // USES assert()

#include <iostream> // TEMPORARY
// ----------------------------------------------------------------------
// These must match the sizes of the buffers in the header file.
const int ALE::CellRefinerTet4::_edgesSize = 9;
const int ALE::CellRefinerTet4::_cellsSize = 36;

// ----------------------------------------------------------------------
// Constructor
ALE::CellRefinerTet4::CellRefinerTet4(const mesh_type& mesh) :
  RefineEdges2(mesh)
{ // constructor
  assert(3 == mesh.getDimension());
} // constructor

// ----------------------------------------------------------------------
// Destructor
ALE::CellRefinerTet4::~CellRefinerTet4(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get number of refined cells for each original cell.
int
ALE::CellRefinerTet4::numNewCells(const point_type cell)
{ // numNewCells
  switch (_cellType(cell)) {
  case TETRAHEDRON:
    return 8;
  case TRIANGLE_COHESIVE_LAGRANGE:
    return 4;
  default:
    assert(0);
    throw ALE::Exception("Unknown cell type.");
  } // switch
} // numNewCells

// ----------------------------------------------------------------------
// Split cell into smaller cells of same type.
void
ALE::CellRefinerTet4::splitCell(const point_type cell,
				const point_type cone[],
				const int coneSize,
				point_type* curNewVertex)
{ // splitCell
  assert(curNewVertex);

  int numEdges = 0;
  const EdgeType* edges;
  
  switch (_cellType(cell)) {
  case TETRAHEDRON:
    _edges_TETRAHEDRON(&edges, &numEdges, cone, coneSize);
    break;
  case TRIANGLE_COHESIVE_LAGRANGE:
    _edges_TRIANGLE_COHESIVE_LAGRANGE(&edges, &numEdges, cone, coneSize);
    break;
  default:
    throw ALE::Exception("Unknown cell type.");
  } // switch

  for(int iEdge=0; iEdge < numEdges; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      // if vertex does not exist
      //std::cout << "Edge: " << edges[iEdge] << ", new vertex: " << *curNewVertex << std::endl;
      _edgeToVertex[edges[iEdge]] = *curNewVertex;
      ++(*curNewVertex);
    } // if
  } // for
} // splitCell

// ----------------------------------------------------------------------
// Split cell into smaller cells of same type.
void
ALE::CellRefinerTet4::splitCellUncensored(const point_type cell,
					  const point_type cone[],
					  const int coneSize,
					  point_type* curNewVertex)
{ // splitCellUncensored
  assert(curNewVertex);

  int numEdges = 0;
  const EdgeType* edges;
  
  const bool uncensored = true;

  switch (_cellType(cell)) {
  case TETRAHEDRON:
    // No censored vertices on normal cell.
    break;
  case TRIANGLE_COHESIVE_LAGRANGE:
    _edges_TRIANGLE_COHESIVE_LAGRANGE(&edges, &numEdges, cone, coneSize, uncensored);
    break;
  default:
    throw ALE::Exception("Unknown cell type.");
  } // switch

  for(int iEdge=0; iEdge < numEdges; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      // if vertex does not exist
      //std::cout << "Edge: " << edges[iEdge] << ", new vertex: " << *curNewVertex << std::endl;
      _edgeToVertex[edges[iEdge]] = *curNewVertex;
      ++(*curNewVertex);
    } // if
  } // for
} // splitCellUncensored

// ----------------------------------------------------------------------
// Get refined cells.
void
ALE::CellRefinerTet4::getNewCells(const point_type** cells,
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
  case TETRAHEDRON: {
    const int coneVertexOffset = orderNewMesh.verticesNormal().min() - orderOldMesh.verticesNormal().min();
    _newCells_TETRAHEDRON(cells, numCells, cone, coneSize, coneVertexOffset);
    break;
  } // TETRAHEDRON
  case TRIANGLE_COHESIVE_LAGRANGE: {
    const int coneVertexOffsetNormal = orderNewMesh.verticesNormal().min() - orderOldMesh.verticesNormal().min();
    const int coneVertexOffsetCensored = orderNewMesh.verticesCensored().min() - orderOldMesh.verticesCensored().min();
    _newCells_TRIANGLE_COHESIVE_LAGRANGE(cells, numCells, cone, coneSize, coneVertexOffsetNormal, coneVertexOffsetCensored);
    break;
  } // TRIANGLE_COHESIVE_LAGRANGE
  default:
    throw ALE::Exception("Unknown cell type.");
  } // switch
} // getNewCells

// ----------------------------------------------------------------------
// Get cell type.
ALE::CellRefinerTet4::CellEnum
ALE::CellRefinerTet4::_cellType(const point_type cell)
{ // _cellType
  assert(!_mesh.getSieve().isNull());

  switch (_mesh.getSieve()->getConeSize(cell)) {
  case 4:
    return TETRAHEDRON;
  case 9:
    return TRIANGLE_COHESIVE_LAGRANGE;
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
ALE::CellRefinerTet4::_edges_TETRAHEDRON(const EdgeType** edges,
					 int* numEdges,
					 const point_type cone[],
					 const int coneSize)
{ // _edges_TETRAHEDRON
  assert(_edgesSize >= 6);
  
  assert(coneSize == 4);
  _edges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
  _edges[1] = EdgeType(std::min(cone[1], cone[2]), std::max(cone[1], cone[2]));
  _edges[2] = EdgeType(std::min(cone[2], cone[0]), std::max(cone[2], cone[0]));
  _edges[3] = EdgeType(std::min(cone[0], cone[3]), std::max(cone[0], cone[3]));
  _edges[4] = EdgeType(std::min(cone[1], cone[3]), std::max(cone[1], cone[3]));
  _edges[5] = EdgeType(std::min(cone[2], cone[3]), std::max(cone[2], cone[3]));
  *numEdges = 6;
  *edges = _edges;
} // _edges_TETRAHEDRON
  
// ----------------------------------------------------------------------
// Get edges of line cohesive cell with Lagrange multipler vertices.
void
ALE::CellRefinerTet4::_edges_TRIANGLE_COHESIVE_LAGRANGE(const EdgeType** edges,
							int* numEdges,
							const point_type cone[],
							const int coneSize,
							const bool uncensored)
{ // _edges_TRIANGLE_COHESIVE_LAGRANGE
  if (uncensored) {
    // Use all vertices
    assert(_edgesSize >= 9);

    assert(coneSize == 9);
    _edges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
    _edges[1] = EdgeType(std::min(cone[1], cone[2]), std::max(cone[1], cone[2]));
    _edges[2] = EdgeType(std::min(cone[2], cone[0]), std::max(cone[2], cone[0]));
    _edges[3] = EdgeType(std::min(cone[3], cone[4]), std::max(cone[3], cone[4]));
    _edges[4] = EdgeType(std::min(cone[4], cone[5]), std::max(cone[4], cone[5]));
    _edges[5] = EdgeType(std::min(cone[5], cone[3]), std::max(cone[5], cone[3]));
    _edges[6] = EdgeType(std::min(cone[6], cone[7]), std::max(cone[6], cone[7]));
    _edges[7] = EdgeType(std::min(cone[7], cone[8]), std::max(cone[7], cone[8]));
    _edges[8] = EdgeType(std::min(cone[8], cone[6]), std::max(cone[8], cone[6]));
    *numEdges = 9;
    *edges = _edges;
  } else {
    // Omit edges with censored (Lagrange multipler) vertices.
    assert(_edgesSize >= 6);

    assert(coneSize == 9);
    _edges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
    _edges[1] = EdgeType(std::min(cone[1], cone[2]), std::max(cone[1], cone[2]));
    _edges[2] = EdgeType(std::min(cone[2], cone[0]), std::max(cone[2], cone[0]));
    _edges[3] = EdgeType(std::min(cone[3], cone[4]), std::max(cone[3], cone[4]));
    _edges[4] = EdgeType(std::min(cone[4], cone[5]), std::max(cone[4], cone[5]));
    _edges[5] = EdgeType(std::min(cone[5], cone[3]), std::max(cone[5], cone[3]));
    *numEdges = 6;
    *edges = _edges;
  } // if/else
} // _edges_TRIANGLE_COHESIVE_LAGRANGE
  
// ----------------------------------------------------------------------
// Get new cells from refinement of a triangular cell.
void
ALE::CellRefinerTet4::_newCells_TETRAHEDRON(const point_type** cells,
					    int *numCells,
					    const point_type cone[],
					    const int coneSize,
					    const int coneVertexOffset)
{ // _newCells_TETRAHEDRON
  const int coneSizeTet4 = 4;
  const int numEdgesTet4 = 6;
  const int numNewCells = 8;
  const int numNewVertices = 6;

  int numEdges = 0;
  const EdgeType* edges;
  _edges_TETRAHEDRON(&edges, &numEdges, cone, coneSize);
  assert(numEdgesTet4 == numEdges);

  assert(_cellsSize >= numNewCells*coneSizeTet4);
  point_type newVertices[numNewVertices];
  for(int iEdge=0, iNewVertex=0; iEdge < numEdgesTet4; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      throw ALE::Exception("Missing edge in refined mesh");
    } // if
    newVertices[iNewVertex++] = _edgeToVertex[edges[iEdge]];
  } // for

  // new cell 0
  _cells[0*4+0] = cone[0]+coneVertexOffset;
  _cells[0*4+1] = newVertices[3];
  _cells[0*4+2] = newVertices[0];
  _cells[0*4+3] = newVertices[2];

  // new cell 1
  _cells[1*4+0] = newVertices[0];
  _cells[1*4+1] = newVertices[1];
  _cells[1*4+2] = newVertices[2];
  _cells[1*4+3] = newVertices[3];

  // new cell 2
  _cells[2*4+0] = newVertices[0];
  _cells[2*4+1] = newVertices[3];
  _cells[2*4+2] = newVertices[4];
  _cells[2*4+3] = newVertices[1];

  // new cell 3
  _cells[3*4+0] = cone[1]+coneVertexOffset;
  _cells[3*4+1] = newVertices[4];
  _cells[3*4+2] = newVertices[1];
  _cells[3*4+3] = newVertices[0];

  // new cell 4
  _cells[4*4+0] = newVertices[2];
  _cells[4*4+1] = newVertices[5];
  _cells[4*4+2] = newVertices[3];
  _cells[4*4+3] = newVertices[1];

  // new cell 5
  _cells[5*4+0] = cone[2]+coneVertexOffset;
  _cells[5*4+1] = newVertices[5];
  _cells[5*4+2] = newVertices[2];
  _cells[5*4+3] = newVertices[1];

  // new cell 6
  _cells[6*4+0] = newVertices[1];
  _cells[6*4+1] = newVertices[4];
  _cells[6*4+2] = newVertices[5];
  _cells[6*4+3] = newVertices[3];

  // new cell 7
  _cells[7*4+0] = cone[3]+coneVertexOffset;
  _cells[7*4+1] = newVertices[3];
  _cells[7*4+2] = newVertices[5];
  _cells[7*4+3] = newVertices[4];

  *numCells = numNewCells;
  *cells = _cells;
} // _newCells_TETRAHEDRON
  
// ----------------------------------------------------------------------
// Get new cells from refinement of a line cohseive cell with Lagrange
// multiplier vertices.
void
ALE::CellRefinerTet4::_newCells_TRIANGLE_COHESIVE_LAGRANGE(const point_type** cells,
							   int *numCells,
							   const point_type cone[],
							   const int coneSize,
							   const int coneVertexOffsetNormal,
							   const int coneVertexOffsetCensored)
{ // _newCells_TRIANGLE_COHESIVE_LAGRANGE
  const int coneSizeTriPrism9 = 9;
  const int numEdgesTriPrism9 = 9;
  const int numNewCells = 4;
  const int numNewVertices = 9;

  int numEdges = 0;
  const EdgeType *edges;
  _edges_TRIANGLE_COHESIVE_LAGRANGE(&edges, &numEdges, cone, coneSize, true);
  assert(numEdgesTriPrism9 == numEdges);

  assert(_cellsSize >= numNewCells*coneSizeTriPrism9);
  point_type newVertices[numNewVertices];
  for(int iEdge=0, iNewVertex=0; iEdge < numEdgesTriPrism9; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      throw ALE::Exception("Missing edge in refined mesh");
    } // if
    newVertices[iNewVertex++] = _edgeToVertex[edges[iEdge]];
  } // for

  _cells[0*9+0] = cone[0]+coneVertexOffsetNormal; // New cell 0
  _cells[0*9+1] = newVertices[0];
  _cells[0*9+2] = newVertices[2];
  _cells[0*9+3] = cone[3]+coneVertexOffsetNormal;
  _cells[0*9+4] = newVertices[3];
  _cells[0*9+5] = newVertices[5];
  _cells[0*9+6] = cone[6]+coneVertexOffsetCensored;
  _cells[0*9+7] = newVertices[6];
  _cells[0*9+8] = newVertices[8];
  
  _cells[1*9+0] = newVertices[0]; // New cell 1
  _cells[1*9+1] = newVertices[1];
  _cells[1*9+2] = newVertices[2];
  _cells[1*9+3] = newVertices[3];
  _cells[1*9+4] = newVertices[4];
  _cells[1*9+5] = newVertices[5];
  _cells[1*9+6] = newVertices[6];
  _cells[1*9+7] = newVertices[7];
  _cells[1*9+8] = newVertices[8];
  
  _cells[2*9+0] = cone[1]+coneVertexOffsetNormal; // New cell 2
  _cells[2*9+1] = newVertices[1];
  _cells[2*9+2] = newVertices[0];
  _cells[2*9+3] = cone[4]+coneVertexOffsetNormal;
  _cells[2*9+4] = newVertices[4];
  _cells[2*9+5] = newVertices[3];
  _cells[2*9+6] = cone[7]+coneVertexOffsetCensored;
  _cells[2*9+7] = newVertices[7];
  _cells[2*9+8] = newVertices[6];
  
  _cells[3*9+0] = cone[2]+coneVertexOffsetNormal; // New cell 3
  _cells[3*9+1] = newVertices[2];
  _cells[3*9+2] = newVertices[1];
  _cells[3*9+3] = cone[5]+coneVertexOffsetNormal;
  _cells[3*9+4] = newVertices[5];
  _cells[3*9+5] = newVertices[4];
  _cells[3*9+6] = cone[8]+coneVertexOffsetCensored;
  _cells[3*9+7] = newVertices[8];
  _cells[3*9+8] = newVertices[7];

  *numCells = numNewCells;
  *cells = _cells;
} // _newCells_TRIANGLE_COHESIVE_LAGRANGE


// End of file 
