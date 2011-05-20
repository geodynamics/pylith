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

#include "CellRefinerHex8.hh" // implementation of class methods

#include "MeshOrder.hh" // USES MeshOrder

#include <cassert> // USES assert()

#include <iostream> // TEMPORARY
// ----------------------------------------------------------------------
// These must match the sizes of the buffers in the header file.
const int ALE::CellRefinerHex8::_edgesSize = 12;
const int ALE::CellRefinerHex8::_facesSize = 6;
const int ALE::CellRefinerHex8::_volumesSize = 1;
const int ALE::CellRefinerHex8::_cellsSize = 64;

// ----------------------------------------------------------------------
// Constructor
ALE::CellRefinerHex8::CellRefinerHex8(const mesh_type& mesh) :
 RefineVol8Face4Edges2(mesh)
{ // constructor
  assert(3 == mesh.getDimension());
} // constructor

// ----------------------------------------------------------------------
// Destructor
ALE::CellRefinerHex8::~CellRefinerHex8(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get number of refined cells for each original cell.
int
ALE::CellRefinerHex8::numNewCells(const point_type cell)
{ // numNewCells
  switch (_cellType(cell)) {
  case HEXAHEDRON:
    return 8;
  case QUAD_COHESIVE_LAGRANGE:
    return 4;
  default:
    assert(0);
    throw ALE::Exception("Unknown cell type.");
  } // switch
} // numNewCells

// ----------------------------------------------------------------------
// Split cell into smaller cells of same type.
void
ALE::CellRefinerHex8::splitCell(const point_type cell,
				const point_type cone[],
				const int coneSize,
				point_type* curNewVertex)
{ // splitCell
  assert(curNewVertex);

  int numEdges = 0;
  const EdgeType* edges;

  int numFaces = 0;
  const FaceType* faces;

  int numVolumes = 0;
  const VolumeType* volumes;

  switch (_cellType(cell)) {
  case HEXAHEDRON:
    _edges_HEXAHEDRON(&edges, &numEdges, cone, coneSize);
    _faces_HEXAHEDRON(&faces, &numFaces, cone, coneSize);
    _volumes_HEXAHEDRON(&volumes, &numVolumes, cone, coneSize);
    break;
  case QUAD_COHESIVE_LAGRANGE:
    _edges_QUAD_COHESIVE_LAGRANGE(&edges, &numEdges, cone, coneSize);
    _faces_QUAD_COHESIVE_LAGRANGE(&faces, &numFaces, cone, coneSize);
    break;
  default:
    assert(0);
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

  for(int iVolume=0; iVolume < numVolumes; ++iVolume) {
    if (_volumeToVertex.find(volumes[iVolume]) == _volumeToVertex.end()) {
      // if vertex does not exist
      //std::cout << "Volume: " << volumes[iVolume] << ", new vertex: " << *curNewVertex << std::endl;
      _volumeToVertex[volumes[iVolume]] = *curNewVertex;
      ++(*curNewVertex);
    } // if
  } // for
} // splitCell

// ----------------------------------------------------------------------
// Split cell into smaller cells of same type.
void
ALE::CellRefinerHex8::splitCellUncensored(const point_type cell,
					  const point_type cone[],
					  const int coneSize,
					  point_type* curNewVertex)
{ // splitCellUncensored
  assert(curNewVertex);

  int numEdges = 0;
  const EdgeType* edges;
  
  int numFaces = 0;
  const FaceType* faces;
  
  const bool uncensored = true;

  switch (_cellType(cell)) {
  case HEXAHEDRON:
    // No censored vertices on normal cell.
    break;
  case QUAD_COHESIVE_LAGRANGE:
    _edges_QUAD_COHESIVE_LAGRANGE(&edges, &numEdges, cone, coneSize, uncensored);
    _faces_QUAD_COHESIVE_LAGRANGE(&faces, &numFaces, cone, coneSize, uncensored);
    break;
  default:
    assert(0);
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

} // splitCellUncensored

// ----------------------------------------------------------------------
// Get refined cells.
void
ALE::CellRefinerHex8::getNewCells(const point_type** cells,
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
  case HEXAHEDRON: {
    const int coneVertexOffset = orderNewMesh.verticesNormal().min() - orderOldMesh.verticesNormal().min();
    _newCells_HEXAHEDRON(cells, numCells, cone, coneSize, coneVertexOffset);
    break;
  } // HEXAHEDRON
  case QUAD_COHESIVE_LAGRANGE: {
    const int coneVertexOffsetNormal = orderNewMesh.verticesNormal().min() - orderOldMesh.verticesNormal().min();
    const int coneVertexOffsetCensored = orderNewMesh.verticesCensored().min() - orderOldMesh.verticesCensored().min();
    _newCells_QUAD_COHESIVE_LAGRANGE(cells, numCells, cone, coneSize, coneVertexOffsetNormal, coneVertexOffsetCensored);
    break;
  } // QUAD_COHESIVE_LAGRANGE
  default:
    assert(0);
    throw ALE::Exception("Unknown cell type.");
  } // switch
} // getNewCells

// ----------------------------------------------------------------------
// Get cell type.
ALE::CellRefinerHex8::CellEnum
ALE::CellRefinerHex8::_cellType(const point_type cell)
{ // _cellType
  assert(!_mesh.getSieve().isNull());

  switch (_mesh.getSieve()->getConeSize(cell)) {
  case 8:
    return HEXAHEDRON;
  case 12:
    return QUAD_COHESIVE_LAGRANGE;
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
// Get edges of hexahedral cell.
void
ALE::CellRefinerHex8::_edges_HEXAHEDRON(const EdgeType** edges,
					int* numEdges,
					const point_type cone[],
					const int coneSize)
{ // _edges_HEXAHEDRON
  assert(_edgesSize >= 12);

  assert(coneSize == 8);
  _edges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
  _edges[1] = EdgeType(std::min(cone[1], cone[2]), std::max(cone[1], cone[2]));
  _edges[2] = EdgeType(std::min(cone[2], cone[3]), std::max(cone[2], cone[3]));
  _edges[3] = EdgeType(std::min(cone[3], cone[0]), std::max(cone[3], cone[0]));

  _edges[4] = EdgeType(std::min(cone[4], cone[5]), std::max(cone[4], cone[5]));
  _edges[5] = EdgeType(std::min(cone[5], cone[6]), std::max(cone[5], cone[6]));
  _edges[6] = EdgeType(std::min(cone[6], cone[7]), std::max(cone[6], cone[7]));
  _edges[7] = EdgeType(std::min(cone[7], cone[4]), std::max(cone[7], cone[4]));

  _edges[8] = EdgeType(std::min(cone[0], cone[4]), std::max(cone[0], cone[4]));
  _edges[9] = EdgeType(std::min(cone[1], cone[5]), std::max(cone[1], cone[5]));
  _edges[10] = EdgeType(std::min(cone[2], cone[6]), std::max(cone[2], cone[6]));
  _edges[11] = EdgeType(std::min(cone[3], cone[7]), std::max(cone[3], cone[7]));

  *numEdges = 12;
  *edges = _edges;
} // _edges_HEXAHEDRON
  
// ----------------------------------------------------------------------
// Get edges of quadrilateral cohesive cell with Lagrange multipler
// vertices.
void
ALE::CellRefinerHex8::_edges_QUAD_COHESIVE_LAGRANGE(const EdgeType** edges,
						    int* numEdges,
						    const point_type cone[],
						    const int coneSize,
						    const bool uncensored)
{ // _edges_QUAD_COHESIVE_LAGRANGE
  if (uncensored) {
    // Include all edges
    assert(_edgesSize >= 12);

    assert(coneSize == 12);
    _edges[0] = EdgeType(std::min(cone[0], cone[1]),
			 std::max(cone[0], cone[1]));
    _edges[1] = EdgeType(std::min(cone[1], cone[2]),
			 std::max(cone[1], cone[2]));
    _edges[2] = EdgeType(std::min(cone[2], cone[3]),
			 std::max(cone[2], cone[3]));
    _edges[3] = EdgeType(std::min(cone[3], cone[0]),
			 std::max(cone[3], cone[0]));

    _edges[4] = EdgeType(std::min(cone[4], cone[5]),
			 std::max(cone[4], cone[5]));
    _edges[5] = EdgeType(std::min(cone[5], cone[6]),
			 std::max(cone[5], cone[6]));
    _edges[6] = EdgeType(std::min(cone[6], cone[7]),
			 std::max(cone[6], cone[7]));
    _edges[7] = EdgeType(std::min(cone[7], cone[4]),
			 std::max(cone[7], cone[4]));

    _edges[8] = EdgeType(std::min(cone[8], cone[9]),
			 std::max(cone[8], cone[9]));
    _edges[9] = EdgeType(std::min(cone[9], cone[10]),
			 std::max(cone[9], cone[10]));
    _edges[10] = EdgeType(std::min(cone[10], cone[11]),
			  std::max(cone[10], cone[11]));
    _edges[11] = EdgeType(std::min(cone[11], cone[8]),
			  std::max(cone[11], cone[8]));

    *numEdges = 12;
    *edges = _edges;
  } else {
    // Omit edges with censored (Lagrange multiplier) vertices.
    assert(_edgesSize >= 8);

    assert(coneSize == 12);
    _edges[0] = EdgeType(std::min(cone[0], cone[1]),
			 std::max(cone[0], cone[1]));
    _edges[1] = EdgeType(std::min(cone[1], cone[2]),
			 std::max(cone[1], cone[2]));
    _edges[2] = EdgeType(std::min(cone[2], cone[3]),
			 std::max(cone[2], cone[3]));
    _edges[3] = EdgeType(std::min(cone[3], cone[0]),
			 std::max(cone[3], cone[0]));

    _edges[4] = EdgeType(std::min(cone[4], cone[5]),
			 std::max(cone[4], cone[5]));
    _edges[5] = EdgeType(std::min(cone[5], cone[6]),
			 std::max(cone[5], cone[6]));
    _edges[6] = EdgeType(std::min(cone[6], cone[7]),
			 std::max(cone[6], cone[7]));
    _edges[7] = EdgeType(std::min(cone[7], cone[4]),
			 std::max(cone[7], cone[4]));

    *numEdges = 8;
    *edges = _edges;
  } // if/else
} // _edges_QUAD_COHESIVE_LAGRANGE
  
// ----------------------------------------------------------------------
// Get faces of hexahedral cell.
void
ALE::CellRefinerHex8::_faces_HEXAHEDRON(const FaceType** faces,
					int* numFaces,
					const point_type cone[],
					const int coneSize)
{ // _faces_HEXAHEDRON
  assert(_facesSize >= 6);
  
  assert(coneSize == 8);
  
  // Face 0, Vertices 0, 1, 5, 4
  if (cone[0] < cone[1] &&
      cone[0] < cone[4] &&
      cone[0] < cone[5]) {
    if (cone[1] < cone[4]) {
      _faces[0] = FaceType(cone[0], cone[1], cone[5], cone[4]);
    } else {
      _faces[0] = FaceType(cone[0], cone[4], cone[5], cone[1]);
    } // if/else
    
  } else if (cone[1] < cone[0] &&
	     cone[1] < cone[4] &&
	     cone[1] < cone[5]) {
    if (cone[0] < cone[5]) {
      _faces[0] = FaceType(cone[1], cone[0], cone[4], cone[5]);
    } else {
      _faces[0] = FaceType(cone[1], cone[5], cone[4], cone[0]);
    } // if/else
    
  } else if (cone[4] < cone[0] &&
	     cone[4] < cone[1] &&
	     cone[4] < cone[5]) {
    if (cone[0] < cone[5]) {
      _faces[0] = FaceType(cone[4], cone[0], cone[1], cone[5]);
    } else {
      _faces[0] = FaceType(cone[4], cone[5], cone[1], cone[0]);
    } // if/else
    
  } else if (cone[5] < cone[0] &&
	     cone[5] < cone[1] &&
	     cone[5] < cone[4]) {
    if (cone[1] < cone[4]) {
      _faces[0] = FaceType(cone[5], cone[1], cone[0], cone[4]);
    } else {
      _faces[0] = FaceType(cone[5], cone[4], cone[0], cone[1]);
    } // if/else
  } else {
    assert(0);
    throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
  } // if/else
  
  // Face 1, Vertices 1, 2, 6, 5
  if (cone[1] < cone[2] &&
      cone[1] < cone[5] &&
      cone[1] < cone[6]) {
    if (cone[2] < cone[5]) {
      _faces[1] = FaceType(cone[1], cone[2], cone[6], cone[5]);
    } else {
      _faces[1] = FaceType(cone[1], cone[5], cone[6], cone[2]);
    } // if/else
    
  } else if (cone[2] < cone[1] &&
	     cone[2] < cone[5] &&
	     cone[2] < cone[6]) {
    if (cone[1] < cone[6]) {
      _faces[1] = FaceType(cone[2], cone[1], cone[5], cone[6]);
    } else {
      _faces[1] = FaceType(cone[2], cone[6], cone[5], cone[1]);
    } // if/else
    
  } else if (cone[5] < cone[1] &&
	     cone[5] < cone[2] &&
	     cone[5] < cone[6]) {
    if (cone[1] < cone[6]) {
      _faces[1] = FaceType(cone[5], cone[1], cone[2], cone[6]);
    } else {
      _faces[1] = FaceType(cone[5], cone[6], cone[2], cone[1]);
    } // if/else
    
  } else if (cone[6] < cone[1] &&
	     cone[6] < cone[2] &&
	     cone[6] < cone[5]) {
    if (cone[2] < cone[5]) {
      _faces[1] = FaceType(cone[6], cone[2], cone[1], cone[5]);
    } else {
      _faces[1] = FaceType(cone[6], cone[5], cone[1], cone[2]);
    } // if/else
  } else {
    assert(0);
    throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
  } // if/else
  
  // Face 2, Vertices 2, 3, 7, 6
  if (cone[2] < cone[3] &&
      cone[2] < cone[6] &&
      cone[2] < cone[7]) {
    if (cone[3] < cone[6]) {
      _faces[2] = FaceType(cone[2], cone[3], cone[7], cone[6]);
    } else {
      _faces[2] = FaceType(cone[2], cone[6], cone[7], cone[3]);
    } // if/else
    
  } else if (cone[3] < cone[2] &&
	     cone[3] < cone[6] &&
	     cone[3] < cone[7]) {
    if (cone[2] < cone[7]) {
      _faces[2] = FaceType(cone[3], cone[2], cone[6], cone[7]);
    } else {
      _faces[2] = FaceType(cone[3], cone[7], cone[6], cone[2]);
    } // if/else
    
  } else if (cone[6] < cone[2] && 
	     cone[6] < cone[3] &&
	     cone[6] < cone[7]) {
    if (cone[2] < cone[7]) {
      _faces[2] = FaceType(cone[6], cone[2], cone[3], cone[7]);
    } else {
      _faces[2] = FaceType(cone[6], cone[7], cone[3], cone[2]);
    } // if/else
    
  } else if (cone[7] < cone[2] &&
	     cone[7] < cone[3] &&
	     cone[7] < cone[6]) {
    if (cone[3] < cone[6]) {
      _faces[2] = FaceType(cone[7], cone[3], cone[2], cone[6]);
    } else {
      _faces[2] = FaceType(cone[7], cone[6], cone[2], cone[3]);
    } // if/else
  } else {
    assert(0);
    throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
  } // if/else
  
  // Face 3, Vertices 0, 3, 7, 4
  if (cone[0] < cone[3] &&
      cone[0] < cone[4] &&
      cone[0] < cone[7]) {
    if (cone[3] < cone[4]) {
      _faces[3] = FaceType(cone[0], cone[3], cone[7], cone[4]);
    } else {
      _faces[3] = FaceType(cone[0], cone[4], cone[7], cone[3]);
    } // if/else
    
  } else if (cone[3] < cone[0] &&
	     cone[3] < cone[4] &&
	     cone[3] < cone[7]) {
    if (cone[0] < cone[7]) {
      _faces[3] = FaceType(cone[3], cone[0], cone[4], cone[7]);
    } else {
      _faces[3] = FaceType(cone[3], cone[7], cone[4], cone[0]);
    } // if/else
    
  } else if (cone[4] < cone[0] && 
	     cone[4] < cone[3] &&
	     cone[4] < cone[7]) {
    if (cone[0] < cone[7]) {
      _faces[3] = FaceType(cone[4], cone[0], cone[3], cone[7]);
    } else {
      _faces[3] = FaceType(cone[4], cone[7], cone[3], cone[0]);
    } // if/else
    
  } else if (cone[7] < cone[0] &&
	     cone[7] < cone[3] &&
	     cone[7] < cone[4]) {
    if (cone[3] < cone[4]) {
      _faces[3] = FaceType(cone[7], cone[3], cone[0], cone[4]);
    } else {
      _faces[3] = FaceType(cone[7], cone[4], cone[0], cone[3]);
    } // if/else
  } else {
    assert(0);
    throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
  } // if/else
  
  // Face 4: Vertices 0, 1, 2, 3
  if (cone[0] < cone[1] &&
      cone[0] < cone[2] &&
      cone[0] < cone[3]) {
    if (cone[1] < cone[3]) {
      _faces[4] = FaceType(cone[0], cone[1], cone[2], cone[3]);
    } else {
      _faces[4] = FaceType(cone[0], cone[3], cone[2], cone[1]);
    } // if/else
    
  } else if (cone[1] < cone[0] &&
	     cone[1] < cone[2] &&
	     cone[1] < cone[3]) {
    if (cone[0] < cone[2]) {
      _faces[4] = FaceType(cone[1], cone[0], cone[3], cone[2]);
    } else {
      _faces[4] = FaceType(cone[1], cone[2], cone[3], cone[0]);
    } // if/else
    
  } else if (cone[2] < cone[0] && 
	     cone[2] < cone[1] &&
	     cone[2] < cone[3]) {
    if (cone[1] < cone[3]) {
      _faces[4] = FaceType(cone[2], cone[1], cone[0], cone[3]);
    } else {
      _faces[4] = FaceType(cone[2], cone[3], cone[0], cone[1]);
    } // if/else
    
  } else if (cone[3] < cone[0] &&
	     cone[3] < cone[1] &&
	     cone[3] < cone[2]) {
    if (cone[0] < cone[2]) {
      _faces[4] = FaceType(cone[3], cone[0], cone[1], cone[2]);
    } else {
      _faces[4] = FaceType(cone[3], cone[2], cone[1], cone[0]);
    } // if/else
  } else {
    assert(0);
    throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
  } // if/else
  
  // Face 5: Vertices 4, 5, 6, 7
  if (cone[4] < cone[5] &&
      cone[4] < cone[6] &&
      cone[4] < cone[7]) {
    if (cone[5] < cone[7]) {
      _faces[5] = FaceType(cone[4], cone[5], cone[6], cone[7]);
    } else {
      _faces[5] = FaceType(cone[4], cone[7], cone[6], cone[5]);
    } // if/else
    
  } else if (cone[5] < cone[4] &&
	     cone[5] < cone[6] &&
	     cone[5] < cone[7]) {
    if (cone[4] < cone[6]) {
      _faces[5] = FaceType(cone[5], cone[4], cone[7], cone[6]);
    } else {
      _faces[5] = FaceType(cone[5], cone[6], cone[7], cone[4]);
    } // if/else
    
  } else if (cone[6] < cone[4] && 
	     cone[6] < cone[5] &&
	     cone[6] < cone[7]) {
    if (cone[5] < cone[7]) {
      _faces[5] = FaceType(cone[6], cone[5], cone[4], cone[7]);
    } else {
      _faces[5] = FaceType(cone[6], cone[7], cone[4], cone[5]);
    } // if/else
    
  } else if (cone[7] < cone[4] &&
	     cone[7] < cone[5] &&
	     cone[7] < cone[6]) {
    if (cone[4] < cone[6]) {
      _faces[5] = FaceType(cone[7], cone[4], cone[5], cone[6]);
    } else {
      _faces[5] = FaceType(cone[7], cone[6], cone[5], cone[4]);
    } // if/else
  } else {
    assert(0);
    throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
  } // if/else
  
  *numFaces = 6;
  *faces = _faces;
} // _faces_HEXAHEDRON
  
// ----------------------------------------------------------------------
// Get faces of quadrilateral cohesive cell with Lagrange multipler
// vertices.
void
ALE::CellRefinerHex8::_faces_QUAD_COHESIVE_LAGRANGE(const FaceType** faces,
						    int* numFaces,
						    const point_type cone[],
						    const int coneSize,
						    const bool uncensored)
{ // _faces_QUAD_COHESIVE_LAGRANGE
  if (uncensored) {
    // Include all faces
    assert(_facesSize >= 3);

    assert(coneSize == 12);

    // Face 0, Vertices 0, 1, 2, 3
    if (cone[0] < cone[1] &&
	cone[0] < cone[2] &&
	cone[0] < cone[3]) {
      if (cone[1] < cone[3]) {
	_faces[0] = FaceType(cone[0], cone[1], cone[2], cone[3]);
      } else {
	_faces[0] = FaceType(cone[0], cone[3], cone[2], cone[1]);
      } // if/else
    
    } else if (cone[1] < cone[0] &&
	       cone[1] < cone[2] &&
	       cone[1] < cone[3]) {
      if (cone[0] < cone[2]) {
	_faces[0] = FaceType(cone[1], cone[0], cone[3], cone[2]);
      } else {
	_faces[0] = FaceType(cone[1], cone[2], cone[3], cone[0]);
      } // if/else
    
    } else if (cone[2] < cone[0] &&
	       cone[2] < cone[1] &&
	       cone[2] < cone[3]) {
      if (cone[1] < cone[3]) {
	_faces[0] = FaceType(cone[2], cone[1], cone[0], cone[3]);
      } else {
	_faces[0] = FaceType(cone[2], cone[3], cone[0], cone[1]);
      } // if/else
    
    } else if (cone[3] < cone[0] &&
	       cone[3] < cone[1] &&
	       cone[3] < cone[2]) {
      if (cone[0] < cone[2]) {
	_faces[0] = FaceType(cone[3], cone[0], cone[1], cone[2]);
      } else {
	_faces[0] = FaceType(cone[3], cone[2], cone[1], cone[0]);
      } // if/else
    } else {
      assert(0);
      throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
    } // if/else
    
    // Face 1: Vertices 4, 5, 6, 7
    if (cone[4] < cone[5] &&
	cone[4] < cone[6] &&
	cone[4] < cone[7]) {
      if (cone[5] < cone[7]) {
	_faces[1] = FaceType(cone[4], cone[5], cone[6], cone[7]);
      } else {
	_faces[1] = FaceType(cone[4], cone[7], cone[6], cone[5]);
      } // if/else
    
    } else if (cone[5] < cone[4] &&
	       cone[5] < cone[6] &&
	       cone[5] < cone[7]) {
      if (cone[4] < cone[6]) {
	_faces[1] = FaceType(cone[5], cone[4], cone[7], cone[6]);
      } else {
	_faces[1] = FaceType(cone[5], cone[6], cone[7], cone[4]);
      } // if/else
    
    } else if (cone[6] < cone[4] && 
	       cone[6] < cone[5] &&
	       cone[6] < cone[7]) {
      if (cone[5] < cone[7]) {
	_faces[1] = FaceType(cone[6], cone[5], cone[4], cone[7]);
      } else {
	_faces[1] = FaceType(cone[6], cone[7], cone[4], cone[5]);
      } // if/else
    
    } else if (cone[7] < cone[4] &&
	       cone[7] < cone[5] &&
	       cone[7] < cone[6]) {
      if (cone[4] < cone[6]) {
	_faces[1] = FaceType(cone[7], cone[4], cone[5], cone[6]);
      } else {
	_faces[1] = FaceType(cone[7], cone[6], cone[5], cone[4]);
      } // if/else
    } else {
      assert(0);
      throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
    } // if/else
  
    // Face 2: Vertices 8, 9, 10, 11
    if (cone[8] < cone[9] &&
	cone[8] < cone[10] &&
	cone[8] < cone[11]) {
      if (cone[9] < cone[11]) {
	_faces[2] = FaceType(cone[8], cone[9], cone[10], cone[11]);
      } else {
	_faces[2] = FaceType(cone[8], cone[11], cone[10], cone[9]);
      } // if/else
    
    } else if (cone[9] < cone[8] &&
	       cone[9] < cone[10] &&
	       cone[9] < cone[11]) {
      if (cone[8] < cone[10]) {
	_faces[2] = FaceType(cone[9], cone[8], cone[11], cone[10]);
      } else {
	_faces[2] = FaceType(cone[9], cone[10], cone[11], cone[8]);
      } // if/else
    
    } else if (cone[10] < cone[8] && 
	       cone[10] < cone[9] &&
	       cone[10] < cone[11]) {
      if (cone[9] < cone[11]) {
	_faces[2] = FaceType(cone[10], cone[9], cone[8], cone[11]);
      } else {
	_faces[2] = FaceType(cone[10], cone[11], cone[8], cone[9]);
      } // if/else
    
    } else if (cone[11] < cone[8] &&
	       cone[11] < cone[9] &&
	       cone[11] < cone[10]) {
      if (cone[8] < cone[10]) {
	_faces[2] = FaceType(cone[11], cone[8], cone[9], cone[10]);
      } else {
	_faces[2] = FaceType(cone[11], cone[10], cone[9], cone[8]);
      } // if/else
    } else {
      assert(0);
      throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
    } // if/else
  
    *numFaces = 3;
    *faces = _faces;
  } else {
    // Omit faces with censored (Lagrange multiplier) vertices.
    assert(_facesSize >= 2);

    assert(coneSize == 12);

    // Face 0, Vertices 0, 1, 2, 3
    if (cone[0] < cone[1] &&
	cone[0] < cone[2] &&
	cone[0] < cone[3]) {
      if (cone[1] < cone[3]) {
	_faces[0] = FaceType(cone[0], cone[1], cone[2], cone[3]);
      } else {
	_faces[0] = FaceType(cone[0], cone[3], cone[2], cone[1]);
      } // if/else
    
    } else if (cone[1] < cone[0] &&
	       cone[1] < cone[2] &&
	       cone[1] < cone[3]) {
      if (cone[0] < cone[2]) {
	_faces[0] = FaceType(cone[1], cone[0], cone[3], cone[2]);
      } else {
	_faces[0] = FaceType(cone[1], cone[2], cone[3], cone[0]);
      } // if/else
    
    } else if (cone[2] < cone[0] &&
	       cone[2] < cone[1] &&
	       cone[2] < cone[3]) {
      if (cone[1] < cone[3]) {
	_faces[0] = FaceType(cone[2], cone[1], cone[0], cone[3]);
      } else {
	_faces[0] = FaceType(cone[2], cone[3], cone[0], cone[1]);
      } // if/else
    
    } else if (cone[3] < cone[0] &&
	       cone[3] < cone[1] &&
	       cone[3] < cone[2]) {
      if (cone[0] < cone[2]) {
	_faces[0] = FaceType(cone[3], cone[0], cone[1], cone[2]);
      } else {
	_faces[0] = FaceType(cone[3], cone[2], cone[1], cone[0]);
      } // if/else
    } else {
      assert(0);
      throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
    } // if/else
    
    // Face 1: Vertices 4, 5, 6, 7
    if (cone[4] < cone[5] &&
	cone[4] < cone[6] &&
	cone[4] < cone[7]) {
      if (cone[5] < cone[7]) {
	_faces[1] = FaceType(cone[4], cone[5], cone[6], cone[7]);
      } else {
	_faces[1] = FaceType(cone[4], cone[7], cone[6], cone[5]);
      } // if/else
    
    } else if (cone[5] < cone[4] &&
	       cone[5] < cone[6] &&
	       cone[5] < cone[7]) {
      if (cone[4] < cone[6]) {
	_faces[1] = FaceType(cone[5], cone[4], cone[7], cone[6]);
      } else {
	_faces[1] = FaceType(cone[5], cone[6], cone[7], cone[4]);
      } // if/else
    
    } else if (cone[6] < cone[4] && 
	       cone[6] < cone[5] &&
	       cone[6] < cone[7]) {
      if (cone[5] < cone[7]) {
	_faces[1] = FaceType(cone[6], cone[5], cone[4], cone[7]);
      } else {
	_faces[1] = FaceType(cone[6], cone[7], cone[4], cone[5]);
      } // if/else
    
    } else if (cone[7] < cone[4] &&
	       cone[7] < cone[5] &&
	       cone[7] < cone[6]) {
      if (cone[4] < cone[6]) {
	_faces[1] = FaceType(cone[7], cone[4], cone[5], cone[6]);
      } else {
	_faces[1] = FaceType(cone[7], cone[6], cone[5], cone[4]);
      } // if/else
    } else {
      assert(0);
      throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
    } // if/else
  
    *numFaces = 2;
    *faces = _faces;
  } // if/else
} // _faces_QUAD_COHESIVE_LAGRANGE
  
// ----------------------------------------------------------------------
// Get volumes of hexahedral cell.
void
ALE::CellRefinerHex8::_volumes_HEXAHEDRON(const VolumeType** volumes,
					int* numVolumes,
					const point_type cone[],
					const int coneSize)
{ // _volumes_HEXAHEDRON
  assert(_volumesSize >= 1);
  
  assert(coneSize == 8);
  
  int sortedCone[8];
  for (int i=0; i < 8; ++i)
    sortedCone[i] = cone[i];
  std::sort(sortedCone, sortedCone+coneSize);

  if (sortedCone[0] == cone[0]) {
    if (cone[1] < cone[3] &&
	cone[1] < cone[4]) {
      _volumes[0] = VolumeType(cone[0], cone[1], cone[2], cone[3],
			       cone[4], cone[5], cone[6], cone[7]);
    } else if (cone[3] < cone[4]) {
      _volumes[0] = VolumeType(cone[0], cone[3], cone[7], cone[4],
			       cone[1], cone[2], cone[6], cone[5]);
    } else {
      _volumes[0] = VolumeType(cone[0], cone[4], cone[5], cone[1],
			       cone[3], cone[7], cone[6], cone[2]);
    } // if/else

  } else if (sortedCone[0] == cone[1]) {
    if (cone[0] < cone[2] &&
	cone[0] < cone[5]) {
      _volumes[0] = VolumeType(cone[1], cone[0], cone[4], cone[5],
			       cone[2], cone[3], cone[7], cone[6]);
    } else if (cone[2] < cone[5]) {
      _volumes[0] = VolumeType(cone[1], cone[2], cone[3], cone[0],
			       cone[5], cone[6], cone[7], cone[4]);
    } else {
      _volumes[0] = VolumeType(cone[1], cone[5], cone[6], cone[2],
			       cone[0], cone[4], cone[7], cone[3]);
    } // if/else

  } else if (sortedCone[0] == cone[2]) {
    if (cone[1] < cone[3] &&
	cone[1] < cone[6]) {
      _volumes[0] = VolumeType(cone[2], cone[1], cone[5], cone[6],
			       cone[3], cone[0], cone[4], cone[7]);
    } else if (cone[3] < cone[6]) {
      _volumes[0] = VolumeType(cone[2], cone[3], cone[0], cone[1],
			       cone[6], cone[7], cone[4], cone[5]);
    } else {
      _volumes[0] = VolumeType(cone[2], cone[6], cone[7], cone[3],
			       cone[1], cone[5], cone[4], cone[0]);
    } // if/else

  } else if (sortedCone[0] == cone[3]) {
    if (cone[0] < cone[2] &&
	cone[0] < cone[7]) {
      _volumes[0] = VolumeType(cone[3], cone[0], cone[1], cone[2],
			       cone[7], cone[4], cone[5], cone[6]);
    } else if (cone[2] < cone[7]) {
      _volumes[0] = VolumeType(cone[3], cone[2], cone[6], cone[7],
			       cone[0], cone[1], cone[5], cone[4]);
    } else {
      _volumes[0] = VolumeType(cone[3], cone[7], cone[4], cone[0],
			       cone[2], cone[6], cone[5], cone[1]);
    } // if/else

  } else if (sortedCone[0] == cone[4]) {
    if (cone[0] < cone[5] &&
	cone[0] < cone[7]) {
      _volumes[0] = VolumeType(cone[4], cone[0], cone[3], cone[7],
			       cone[5], cone[1], cone[2], cone[6]);
    } else if (cone[5] < cone[7]) {
      _volumes[0] = VolumeType(cone[4], cone[5], cone[1], cone[0],
			       cone[7], cone[6], cone[2], cone[3]);
    } else {
      _volumes[0] = VolumeType(cone[4], cone[7], cone[6], cone[5],
			       cone[0], cone[3], cone[2], cone[1]);
    } // if/else

  } else if (sortedCone[0] == cone[5]) {
    if (cone[1] < cone[4] &&
	cone[1] < cone[6]) {
      _volumes[0] = VolumeType(cone[5], cone[1], cone[0], cone[4],
			       cone[6], cone[2], cone[3], cone[7]);
    } else if (cone[4] < cone[6]) {
      _volumes[0] = VolumeType(cone[5], cone[4], cone[7], cone[6],
			       cone[1], cone[0], cone[3], cone[2]);
    } else {
      _volumes[0] = VolumeType(cone[5], cone[6], cone[2], cone[1],
			       cone[4], cone[7], cone[3], cone[0]);
    } // if/else


  } else if (sortedCone[0] == cone[6]) {
    if (cone[2] < cone[5] &&
	cone[2] < cone[7]) {
      _volumes[0] = VolumeType(cone[6], cone[2], cone[1], cone[5],
			       cone[7], cone[3], cone[0], cone[4]);
    } else if (cone[5] < cone[7]) {
      _volumes[0] = VolumeType(cone[6], cone[5], cone[4], cone[7],
			       cone[2], cone[1], cone[0], cone[3]);
    } else {
      _volumes[0] = VolumeType(cone[6], cone[7], cone[3], cone[2],
			       cone[5], cone[4], cone[0], cone[1]);
    } // if/else


  } else if (sortedCone[0] == cone[7]) {
    if (cone[3] < cone[4] &&
	cone[3] < cone[6]) {
      _volumes[0] = VolumeType(cone[7], cone[3], cone[2], cone[6],
			       cone[4], cone[0], cone[1], cone[5]);
    } else if (cone[4] < cone[6]) {
      _volumes[0] = VolumeType(cone[7], cone[4], cone[0], cone[3],
			       cone[6], cone[5], cone[1], cone[2]);
    } else {
      _volumes[0] = VolumeType(cone[7], cone[6], cone[5], cone[4],
			       cone[3], cone[2], cone[1], cone[0]);
    } // if/else

  } else {
    assert(0);
    throw ALE::Exception("Could not determine quad face orientation during uniform global refinement.");
  } // if/else
  *numVolumes = 1;
  *volumes = _volumes;
} // _volumes_HEXAHEDRON
  
// ----------------------------------------------------------------------
// Get new cells from refinement of a hexahedral cell.
void
ALE::CellRefinerHex8::_newCells_HEXAHEDRON(const point_type** cells,
					   int *numCells,
					   const point_type cone[],
					   const int coneSize,
					   const int coneVertexOffset)
{ // _newCells_HEXAHEDRON
  const int coneSizeHex8 = 8;
  const int numEdgesHex8 = 12;
  const int numFacesHex8 = 6;
  const int numVolumesHex8 = 1;
  const int numNewCells = 8;
  const int numNewVertices = 19;

  int numEdges = 0;
  const EdgeType* edges;
  _edges_HEXAHEDRON(&edges, &numEdges, cone, coneSize);
  assert(numEdgesHex8 == numEdges);

  int numFaces = 0;
  const FaceType* faces;
  _faces_HEXAHEDRON(&faces, &numFaces, cone, coneSize);
  assert(numFacesHex8 == numFaces);

  int numVolumes = 0;
  const VolumeType* volumes;
  _volumes_HEXAHEDRON(&volumes, &numVolumes, cone, coneSize);
  assert(numVolumesHex8 == numVolumes);

  assert(_cellsSize >= numNewCells*coneSizeHex8);
  point_type newVertices[numNewVertices];
  int iNewVertex = 0;
  for(int iEdge=0; iEdge < numEdgesHex8; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      throw ALE::Exception("Missing edge in refined mesh");
    } // if
    newVertices[iNewVertex++] = _edgeToVertex[edges[iEdge]];
  } // for
  for(int iFace=0; iFace < numFacesHex8; ++iFace) {
    if (_faceToVertex.find(faces[iFace]) == _faceToVertex.end()) {
      throw ALE::Exception("Missing face in refined mesh");
    } // if
    newVertices[iNewVertex++] = _faceToVertex[faces[iFace]];
  } // for
  for(int iVolume=0; iVolume < numVolumesHex8; ++iVolume) {
    if (_volumeToVertex.find(volumes[iVolume]) == _volumeToVertex.end()) {
      throw ALE::Exception("Missing volume in refined mesh");
    } // if
    newVertices[iNewVertex++] = _volumeToVertex[volumes[iVolume]];
  } // for

  // new cell 0
  _cells[0*8+0] = cone[0] + coneVertexOffset;
  _cells[0*8+1] = newVertices[ 0];
  _cells[0*8+2] = newVertices[16];
  _cells[0*8+3] = newVertices[ 3];
  _cells[0*8+4] = newVertices[ 8];
  _cells[0*8+5] = newVertices[12];
  _cells[0*8+6] = newVertices[18];
  _cells[0*8+7] = newVertices[15];

  // new cell 1
  _cells[1*8+0] = cone[1] + coneVertexOffset;
  _cells[1*8+1] = newVertices[ 1];
  _cells[1*8+2] = newVertices[16];
  _cells[1*8+3] = newVertices[ 0];
  _cells[1*8+4] = newVertices[ 9];
  _cells[1*8+5] = newVertices[13];
  _cells[1*8+6] = newVertices[18];
  _cells[1*8+7] = newVertices[12];

  // new cell 2
  _cells[2*8+0] = cone[2] + coneVertexOffset;
  _cells[2*8+1] = newVertices[ 2];
  _cells[2*8+2] = newVertices[16];
  _cells[2*8+3] = newVertices[ 1];
  _cells[2*8+4] = newVertices[10];
  _cells[2*8+5] = newVertices[14];
  _cells[2*8+6] = newVertices[18];
  _cells[2*8+7] = newVertices[13];

  // new cell 3
  _cells[3*8+0] = cone[3] + coneVertexOffset;
  _cells[3*8+1] = newVertices[ 3];
  _cells[3*8+2] = newVertices[16];
  _cells[3*8+3] = newVertices[ 2];
  _cells[3*8+4] = newVertices[11];
  _cells[3*8+5] = newVertices[15];
  _cells[3*8+6] = newVertices[18];
  _cells[3*8+7] = newVertices[14];

  // new cell 4
  _cells[4*8+0] = newVertices[ 8];
  _cells[4*8+1] = newVertices[12];
  _cells[4*8+2] = newVertices[18];
  _cells[4*8+3] = newVertices[15];
  _cells[4*8+4] = cone[4] + coneVertexOffset;
  _cells[4*8+5] = newVertices[ 4];
  _cells[4*8+6] = newVertices[17];
  _cells[4*8+7] = newVertices[ 7];

  // new cell 5
  _cells[5*8+0] = newVertices[ 9];
  _cells[5*8+1] = newVertices[13];
  _cells[5*8+2] = newVertices[18];
  _cells[5*8+3] = newVertices[12];
  _cells[5*8+4] = cone[5] + coneVertexOffset;
  _cells[5*8+5] = newVertices[ 5];
  _cells[5*8+6] = newVertices[17];
  _cells[5*8+7] = newVertices[ 4];

  // new cell 6
  _cells[6*8+0] = newVertices[10];
  _cells[6*8+1] = newVertices[14];
  _cells[6*8+2] = newVertices[18];
  _cells[6*8+3] = newVertices[13];
  _cells[6*8+4] = cone[6] + coneVertexOffset;
  _cells[6*8+5] = newVertices[ 6];
  _cells[6*8+6] = newVertices[17];
  _cells[6*8+7] = newVertices[ 5];

  // new cell 7
  _cells[7*8+0] = newVertices[11];
  _cells[7*8+1] = newVertices[15];
  _cells[7*8+2] = newVertices[18];
  _cells[7*8+3] = newVertices[14];
  _cells[7*8+4] = cone[7] + coneVertexOffset;
  _cells[7*8+5] = newVertices[ 7];
  _cells[7*8+6] = newVertices[17];
  _cells[7*8+7] = newVertices[ 6];

  *numCells = numNewCells;
  *cells = _cells;
} // _newCells_QUADRILATERAL
  
// ----------------------------------------------------------------------
// Get new cells from refinement of a quad cohseive cell with Lagrange
// multiplier vertices.
void
ALE::CellRefinerHex8::_newCells_QUAD_COHESIVE_LAGRANGE(const point_type** cells,
						       int *numCells,
						       const point_type cone[],
						       const int coneSize,
						       const int coneVertexOffsetNormal,
						       const int coneVertexOffsetCensored)
{ // _newCells_QUAD_COHESIVE_LAGRANGE
  const int coneSizeQuad12 = 12;
  const int numEdgesQuad12 = 12;
  const int numFacesQuad12 = 3;
  const int numNewCells = 4;
  const int numNewVertices = 15;

  int numEdges = 0;
  const EdgeType *edges;
  _edges_QUAD_COHESIVE_LAGRANGE(&edges, &numEdges, cone, coneSize, true);
  assert(numEdgesQuad12 == numEdges);

  int numFaces = 0;
  const FaceType *faces;
  _faces_QUAD_COHESIVE_LAGRANGE(&faces, &numFaces, cone, coneSize, true);
  assert(numFacesQuad12 == numFaces);

  assert(_cellsSize >= numNewCells*coneSizeQuad12);
  point_type newVertices[numNewVertices];
  int iNewVertex = 0;
  for(int iEdge=0; iEdge < numEdgesQuad12; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      throw ALE::Exception("Missing edge in refined mesh");
    } // if
    newVertices[iNewVertex++] = _edgeToVertex[edges[iEdge]];
  } // for
  for(int iFace=0; iFace < numFacesQuad12; ++iFace) {
    if (_faceToVertex.find(faces[iFace]) == _faceToVertex.end()) {
      throw ALE::Exception("Missing face in refined mesh");
    } // if
    newVertices[iNewVertex++] = _faceToVertex[faces[iFace]];
  } // for
  assert(numNewVertices == iNewVertex);

  // new cell 0
  _cells[0*12+ 0] = cone[0] + coneVertexOffsetNormal;
  _cells[0*12+ 1] = newVertices[ 0];
  _cells[0*12+ 2] = newVertices[12];
  _cells[0*12+ 3] = newVertices[ 3];
  _cells[0*12+ 4] = cone[4] + coneVertexOffsetNormal;
  _cells[0*12+ 5] = newVertices[ 4];
  _cells[0*12+ 6] = newVertices[13];
  _cells[0*12+ 7] = newVertices[ 7];
  _cells[0*12+ 8] = cone[8] + coneVertexOffsetCensored;
  _cells[0*12+ 9] = newVertices[ 8];
  _cells[0*12+10] = newVertices[14];
  _cells[0*12+11] = newVertices[11];

  // new cell 1
  _cells[1*12+ 0] = cone[1] + coneVertexOffsetNormal;
  _cells[1*12+ 1] = newVertices[ 1];
  _cells[1*12+ 2] = newVertices[12];
  _cells[1*12+ 3] = newVertices[ 0];
  _cells[1*12+ 4] = cone[5] + coneVertexOffsetNormal;
  _cells[1*12+ 5] = newVertices[ 5];
  _cells[1*12+ 6] = newVertices[13];
  _cells[1*12+ 7] = newVertices[ 4];
  _cells[1*12+ 8] = cone[9] + coneVertexOffsetCensored;
  _cells[1*12+ 9] = newVertices[ 9];
  _cells[1*12+10] = newVertices[14];
  _cells[1*12+11] = newVertices[ 8];

  // new cell 2
  _cells[2*12+ 0] = cone[2] + coneVertexOffsetNormal;
  _cells[2*12+ 1] = newVertices[ 2];
  _cells[2*12+ 2] = newVertices[12];
  _cells[2*12+ 3] = newVertices[ 1];
  _cells[2*12+ 4] = cone[6] + coneVertexOffsetNormal;
  _cells[2*12+ 5] = newVertices[ 6];
  _cells[2*12+ 6] = newVertices[13];
  _cells[2*12+ 7] = newVertices[ 5];
  _cells[2*12+ 8] = cone[10] + coneVertexOffsetCensored;
  _cells[2*12+ 9] = newVertices[10];
  _cells[2*12+10] = newVertices[14];
  _cells[2*12+11] = newVertices[ 9];

  // new cell 3
  _cells[3*12+ 0] = cone[3] + coneVertexOffsetNormal;
  _cells[3*12+ 1] = newVertices[ 3];
  _cells[3*12+ 2] = newVertices[12];
  _cells[3*12+ 3] = newVertices[ 2];
  _cells[3*12+ 4] = cone[7] + coneVertexOffsetNormal;
  _cells[3*12+ 5] = newVertices[ 7];
  _cells[3*12+ 6] = newVertices[13];
  _cells[3*12+ 7] = newVertices[ 6];
  _cells[3*12+ 8] = cone[11] + coneVertexOffsetCensored;
  _cells[3*12+ 9] = newVertices[11];
  _cells[3*12+10] = newVertices[14];
  _cells[3*12+11] = newVertices[10];

  *numCells = 4;
  *cells = _cells;
} // _newCells_QUAD_COHESIVE_LAGRANGE


// End of file 
