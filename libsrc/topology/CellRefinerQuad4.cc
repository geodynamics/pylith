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

#include "CellRefinerQuad4.hh" // implementation of class methods

#include "MeshOrder.hh" // USES MeshOrder

#include <cassert> // USES assert()

#include <iostream> // TEMPORARY
// ----------------------------------------------------------------------
// These must match the sizes of the buffers in the header file.
const int ALE::CellRefinerQuad4::_edgesSize = 4;
const int ALE::CellRefinerQuad4::_facesSize = 1;
const int ALE::CellRefinerQuad4::_cellsSize = 16;

// ----------------------------------------------------------------------
// Constructor
ALE::CellRefinerQuad4::CellRefinerQuad4(const mesh_type& mesh) :
 RefineFace4Edges2(mesh)
{ // constructor
  assert(2 == mesh.getDimension());
} // constructor

// ----------------------------------------------------------------------
// Destructor
ALE::CellRefinerQuad4::~CellRefinerQuad4(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get number of refined cells for each original cell.
int
ALE::CellRefinerQuad4::numNewCells(const point_type cell)
{ // numNewCells
  switch (_cellType(cell)) {
  case QUADRILATERAL:
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
ALE::CellRefinerQuad4::splitCell(const point_type cell,
				const point_type cone[],
				const int coneSize,
				point_type* curNewVertex)
{ // splitCell
  assert(curNewVertex);

  int numEdges = 0;
  const EdgeType* edges;

  int numFaces = 0;
  const FaceType* faces;

  switch (_cellType(cell)) {
  case QUADRILATERAL:
    _edges_QUADRILATERAL(&edges, &numEdges, cone, coneSize);
    _faces_QUADRILATERAL(&faces, &numFaces, cone, coneSize);
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
      //std::cout << "Edge: " << edges[iEdge] << ", new vertex: " << *curNewVertex << std::endl;
      _edgeToVertex[edges[iEdge]] = *curNewVertex;
      ++(*curNewVertex);
    } // if
  } // for

  for(int iFace=0; iFace < numFaces; ++iFace) {
    if (_faceToVertex.find(faces[iFace]) == _faceToVertex.end()) {
      // if vertex does not exist
      //std::cout << "Face: " << faces[iFace] << ", new vertex: " << *curNewVertex << std::endl;
      _faceToVertex[faces[iFace]] = *curNewVertex;
      ++(*curNewVertex);
    } // if
  } // for
} // splitCell

// ----------------------------------------------------------------------
// Split cell into smaller cells of same type.
void
ALE::CellRefinerQuad4::splitCellUncensored(const point_type cell,
					  const point_type cone[],
					  const int coneSize,
					  point_type* curNewVertex)
{ // splitCellUncensored
  assert(curNewVertex);

  int numEdges = 0;
  const EdgeType* edges;
  
  const bool uncensored = true;

  switch (_cellType(cell)) {
  case QUADRILATERAL:
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
      //std::cout << "Edge: " << edges[iEdge] << ", new vertex: " << *curNewVertex << std::endl;
      _edgeToVertex[edges[iEdge]] = *curNewVertex;
      ++(*curNewVertex);
    } // if
  } // for

} // splitCellUncensored

// ----------------------------------------------------------------------
// Get refined cells.
void
ALE::CellRefinerQuad4::getNewCells(const point_type** cells,
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
  case QUADRILATERAL: {
    const int coneVertexOffset = orderNewMesh.verticesNormal().min() - orderOldMesh.verticesNormal().min();
    _newCells_QUADRILATERAL(cells, numCells, cone, coneSize, coneVertexOffset);
    break;
  } // QUADRILATERAL
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
ALE::CellRefinerQuad4::CellEnum
ALE::CellRefinerQuad4::_cellType(const point_type cell)
{ // _cellType
  assert(!_mesh.getSieve().isNull());

  switch (_mesh.getSieve()->getConeSize(cell)) {
  case 4:
    return QUADRILATERAL;
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
// Get edges of quadrilateral cell.
void
ALE::CellRefinerQuad4::_edges_QUADRILATERAL(const EdgeType** edges,
					    int* numEdges,
					    const point_type cone[],
					    const int coneSize)
{ // _edges_QUADRILATERAL
  assert(_edgesSize >= 4);

  assert(coneSize == 4);
  _edges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
  _edges[1] = EdgeType(std::min(cone[1], cone[2]), std::max(cone[1], cone[2]));
  _edges[2] = EdgeType(std::min(cone[2], cone[3]), std::max(cone[2], cone[3]));
  _edges[3] = EdgeType(std::min(cone[3], cone[0]), std::max(cone[3], cone[0]));
  *numEdges = 4;
  *edges = _edges;
} // _edges_QUADRILATERAL
  
// ----------------------------------------------------------------------
// Get edges of line cohesive cell with Lagrange multipler vertices.
void
ALE::CellRefinerQuad4::_edges_LINE_COHESIVE_LAGRANGE(const EdgeType** edges,
						    int* numEdges,
						    const point_type cone[],
						    const int coneSize,
						    const bool uncensored)
{ // _edges_LINE_COHESIVE_LAGRANGE
  if (uncensored) {
    // Include all edges
    assert(_edgesSize >= 3);

    assert(coneSize == 6);
    _edges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
    _edges[1] = EdgeType(std::min(cone[2], cone[3]), std::max(cone[2], cone[3]));
    _edges[2] = EdgeType(std::min(cone[4], cone[5]), std::max(cone[4], cone[5]));
    *numEdges = 3;
    *edges = _edges;
  } else {
    // Omit edges with censored (Lagrange multiplier) vertices.
    assert(_edgesSize >= 2);

    assert(coneSize == 6);
    _edges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
    _edges[1] = EdgeType(std::min(cone[2], cone[3]), std::max(cone[2], cone[3]));
    *numEdges = 2;
    *edges = _edges;
  } // if/else
} // _edges_LINE_COHESIVE_LAGRANGE
  
// ----------------------------------------------------------------------
// Get faces of quadrilateral cell.
void
ALE::CellRefinerQuad4::_faces_QUADRILATERAL(const FaceType** faces,
					    int* numFaces,
					    const point_type cone[],
					    const int coneSize)
{ // _faces_QUADRILATERAL
  assert(_facesSize >= 1);
  
  assert(coneSize == 4);
  
  int sortedCone[4];
  for (int i=0; i < 4; ++i)
    sortedCone[i] = cone[i];
  std::sort(sortedCone, sortedCone+coneSize);

  if (sortedCone[0] == cone[0]) {
    if (cone[1] < cone[3]) {
      _faces[0] = FaceType(cone[0], cone[1], cone[2], cone[3]);
    } else {
      _faces[0] = FaceType(cone[0], cone[3], cone[2], cone[1]);
    } // if/else

  } else if (sortedCone[0] == cone[1]) {
    if (cone[2] < cone[0]) {
      _faces[0] = FaceType(cone[1], cone[2], cone[3], cone[0]);
    } else {
      _faces[0] = FaceType(cone[1], cone[0], cone[3], cone[2]);
    } // if/else

  } else if (sortedCone[0] == cone[2]) {
    if (cone[3] < cone[1]) {
      _faces[0] = FaceType(cone[2], cone[3], cone[0], cone[1]);
    } else {
      _faces[0] = FaceType(cone[2], cone[1], cone[3], cone[0]);
    } // if/else

  } else if (sortedCone[0] == cone[3]) {
    if (cone[0] < cone[2]) {
      _faces[0] = FaceType(cone[3], cone[0], cone[1], cone[2]);
    } else {
      _faces[0] = FaceType(cone[3], cone[2], cone[1], cone[0]);
    } // if/else
  } else {
    assert(0);
    throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
  } // if/else
  *numFaces = 1;
  *faces = _faces;
} // _faces_QUADRILATERAL
  
// ----------------------------------------------------------------------
// Get new cells from refinement of a triangular cell.
void
  ALE::CellRefinerQuad4::_newCells_QUADRILATERAL(const point_type** cells,
						 int *numCells,
						 const point_type cone[],
						 const int coneSize,
						 const int coneVertexOffset)
{ // _newCells_QUADRILATERAL
  const int coneSizeQuad4 = 4;
  const int numEdgesQuad4 = 4;
  const int numFacesQuad4 = 1;
  const int numNewCells = 4;
  const int numNewVertices = 5;

  int numEdges = 0;
  const EdgeType* edges;
  _edges_QUADRILATERAL(&edges, &numEdges, cone, coneSize);
  assert(numEdgesQuad4 == numEdges);

  int numFaces = 0;
  const FaceType* faces;
  _faces_QUADRILATERAL(&faces, &numFaces, cone, coneSize);
  assert(numFacesQuad4 == numFaces);

  assert(_cellsSize >= numNewCells*coneSizeQuad4);
  point_type newVertices[numNewVertices];
  int iNewVertex = 0;
  for(int iEdge=0; iEdge < numEdgesQuad4; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      throw ALE::Exception("Missing edge in refined mesh");
    } // if
    newVertices[iNewVertex++] = _edgeToVertex[edges[iEdge]];
  } // for
  for(int iFace=0; iFace < numFacesQuad4; ++iFace) {
    if (_faceToVertex.find(faces[iFace]) == _faceToVertex.end()) {
      throw ALE::Exception("Missing face in refined mesh");
    } // if
    newVertices[iNewVertex++] = _faceToVertex[faces[iFace]];
  } // for

  // new cell 0
  _cells[0*4+0] = cone[0] + coneVertexOffset;
  _cells[0*4+1] = newVertices[0];
  _cells[0*4+2] = newVertices[4];
  _cells[0*4+3] = newVertices[3];

  // new cell 1
  _cells[1*4+0] = cone[1] + coneVertexOffset;
  _cells[1*4+1] = newVertices[1];
  _cells[1*4+2] = newVertices[4];
  _cells[1*4+3] = newVertices[0];

  // new cell 2
  _cells[2*4+0] = cone[3] + coneVertexOffset;
  _cells[2*4+1] = newVertices[3];
  _cells[2*4+2] = newVertices[4];
  _cells[2*4+3] = newVertices[2];

  // new cell 3
  _cells[3*4+0] = cone[2] + coneVertexOffset;
  _cells[3*4+1] = newVertices[2];
  _cells[3*4+2] = newVertices[4];
  _cells[3*4+3] = newVertices[1];

  *numCells = numNewCells;
  *cells = _cells;
} // _newCells_QUADRILATERAL
  
// ----------------------------------------------------------------------
// Get new cells from refinement of a line cohseive cell with Lagrange
// multiplier vertices.
void
ALE::CellRefinerQuad4::_newCells_LINE_COHESIVE_LAGRANGE(const point_type** cells,
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

  assert(_cellsSize >= numNewCells*coneSizeLine6);
  point_type newVertices[numNewVertices];
  for(int iEdge=0, iNewVertex=0; iEdge < numEdgesLine6; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      throw ALE::Exception("Missing edge in refined mesh");
    } // if
    newVertices[iNewVertex++] = _edgeToVertex[edges[iEdge]];
  } // for

  // new cell 0
  _cells[0*6+0] = cone[0] + coneVertexOffsetNormal;
  _cells[0*6+1] = newVertices[0];
  _cells[0*6+2] = cone[2] + coneVertexOffsetNormal;
  _cells[0*6+3] = newVertices[1];
  _cells[0*6+4] = cone[4] + coneVertexOffsetCensored;
  _cells[0*6+5] = newVertices[2];

  // new cell 1
  _cells[1*6+0] = newVertices[0];
  _cells[1*6+1] = cone[1] + coneVertexOffsetNormal;
  _cells[1*6+2] = newVertices[1];
  _cells[1*6+3] = cone[3] + coneVertexOffsetNormal;
  _cells[1*6+4] = newVertices[2];
  _cells[1*6+5] = cone[5] + coneVertexOffsetCensored;
  
  *numCells = 2;
  *cells = _cells;
} // _newCells_LINE_COHESIVE_LAGRANGE


// End of file 
