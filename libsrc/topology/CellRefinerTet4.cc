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

#include "CellRefinerTet4.hh" // implementation of class methods

#include "MeshOrder.hh" // USES MeshOrder

#include <cassert> // USES assert()

#include <iostream> // TEMPORARY
// ----------------------------------------------------------------------
// Constructor
ALE::CellRefinerTet4::CellRefinerTet4(const mesh_type& mesh) :
  _mesh(mesh)
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
      std::cout << "Edge: " << edges[iEdge] << ", new vertex: " << *curNewVertex << std::endl;
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
      std::cout << "Edge: " << edges[iEdge] << ", new vertex: " << *curNewVertex << std::endl;
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
// Set coordinates of new vertices.
void
ALE::CellRefinerTet4::setCoordsNewVertices(const ALE::Obj<mesh_type::real_section_type>& newCoordsSection,
					   const ALE::Obj<mesh_type::real_section_type>& oldCoordsSection)
{ // setCoordsNewVertices
  assert(!newCoordsSection.isNull());
  assert(!oldCoordsSection.isNull());

  double coordinatesVertex[3];

  assert(_edgeToVertex.size() > 0);
  const int spaceDim = newCoordsSection->getFiberDimension(_edgeToVertex.begin()->second);
  assert(spaceDim > 0 && spaceDim <= 3);

  const edge_map_type::const_iterator edgesEnd = _edgeToVertex.end();
  for (edge_map_type::const_iterator e_iter = _edgeToVertex.begin(); e_iter != edgesEnd; ++e_iter) {
    const point_type newVertex = e_iter->second;
    const point_type edgeVertexA = e_iter->first.first;
    const point_type edgeVertexB = e_iter->first.second;

    assert(spaceDim == oldCoordsSection->getFiberDimension(edgeVertexA));
    assert(spaceDim == oldCoordsSection->getFiberDimension(edgeVertexB));
    assert(spaceDim == newCoordsSection->getFiberDimension(newVertex));

    const mesh_type::real_section_type::value_type* coordsA = oldCoordsSection->restrictPoint(edgeVertexA);
    const mesh_type::real_section_type::value_type* coordsB = oldCoordsSection->restrictPoint(edgeVertexB);
    for (int i=0; i < spaceDim; ++i)
      coordinatesVertex[i] = 0.5*(coordsA[i] + coordsB[i]);

    newCoordsSection->updatePoint(newVertex, coordinatesVertex);
  } // for
} // setCoordsNewVertices

// ----------------------------------------------------------------------
// Add space for new vertices in group.
void
ALE::CellRefinerTet4::groupAddNewVertices(const ALE::Obj<mesh_type::int_section_type>& newGroup,
					  const ALE::Obj<mesh_type::int_section_type>& oldGroup)
{ // groupAddNewVertices
  assert(!newGroup.isNull());
  assert(!oldGroup.isNull());

  const edge_map_type::const_iterator edgesEnd = _edgeToVertex.end();
  for (edge_map_type::const_iterator e_iter = _edgeToVertex.begin(); e_iter != edgesEnd; ++e_iter) {
    const point_type newVertex = e_iter->second;
    const point_type edgeVertexA = e_iter->first.first;
    const point_type edgeVertexB = e_iter->first.second;

    if (oldGroup->getFiberDimension(edgeVertexA) && oldGroup->getFiberDimension(edgeVertexB)) {
      if (oldGroup->restrictPoint(edgeVertexA)[0] == oldGroup->restrictPoint(edgeVertexB)[0]) {
	  newGroup->setFiberDimension(newVertex, 1);
      } // if
    } // if
  } // for
} // groupAddNewVertices

// ----------------------------------------------------------------------
// Set new vertices in group.
void
ALE::CellRefinerTet4::groupSetNewVertices(const ALE::Obj<mesh_type::int_section_type>& newGroup,
					  const ALE::Obj<mesh_type::int_section_type>& oldGroup)
{ // groupSetNewVertices
  assert(!newGroup.isNull());
  assert(!oldGroup.isNull());

  const edge_map_type::const_iterator edgesEnd = _edgeToVertex.end();
  for (edge_map_type::const_iterator e_iter = _edgeToVertex.begin(); e_iter != edgesEnd; ++e_iter) {
    const point_type newVertex = e_iter->second;
    const point_type edgeVertexA = e_iter->first.first;
    const point_type edgeVertexB = e_iter->first.second;

    if (oldGroup->getFiberDimension(edgeVertexA) && oldGroup->getFiberDimension(edgeVertexB)) {
      if (oldGroup->restrictPoint(edgeVertexA)[0] == oldGroup->restrictPoint(edgeVertexB)[0]) {
	newGroup->updatePoint(newVertex, oldGroup->restrictPoint(edgeVertexA));
	std::cout << "Adding new vertex: " << newVertex << " based on old vertices " << edgeVertexA << " and " << edgeVertexB << std::endl;
      } // if
    } // if
  } // for
} // groupSetNewVertices

// ----------------------------------------------------------------------
// Add new vertices to label.
void
ALE::CellRefinerTet4::labelAddNewVertices(const ALE::Obj<mesh_type>& newMesh,
					  const ALE::Obj<mesh_type>& oldMesh,
					  const char* labelName)
{ // labelAddNewVertices
  assert(!newMesh.isNull());
  assert(!oldMesh.isNull());

  const Obj<mesh_type::label_sequence>& oldLabelVertices = oldMesh->getLabelStratum(labelName, 0);
  assert(!oldLabelVertices.isNull());

  const Obj<mesh_type::label_type>& oldLabel = oldMesh->getLabel(labelName);
  assert(!oldLabel.isNull());
  const Obj<mesh_type::label_type>& newLabel = newMesh->getLabel(labelName);
  assert(!newLabel.isNull());

  const int defaultValue = -999;

  const edge_map_type::const_iterator edgesEnd = _edgeToVertex.end();
  for (edge_map_type::const_iterator e_iter = _edgeToVertex.begin(); e_iter != edgesEnd; ++e_iter) {
    const point_type newVertex = e_iter->second;
    const point_type edgeVertexA = e_iter->first.first;
    const point_type edgeVertexB = e_iter->first.second;

    const int valueA = oldMesh->getValue(oldLabel, edgeVertexA, defaultValue);
    const int valueB = oldMesh->getValue(oldLabel, edgeVertexB, defaultValue);

    if (valueA != defaultValue && valueA == valueB) {
      newMesh->setValue(newLabel, newVertex, valueA);
    } // if
  } // for
} // labelAddNewVertices

// ----------------------------------------------------------------------
// Calculate new overlap.
void
ALE::CellRefinerTet4::overlapAddNewVertices(const Obj<mesh_type>& newMesh,
					    const MeshOrder& orderNewMesh,
					    const Obj<mesh_type>& oldMesh,
					    const MeshOrder& orderOldMesh)
{ // overlapAddNewVertices
  assert(!newMesh.isNull());
  assert(!oldMesh.isNull());

  Obj<mesh_type::send_overlap_type> newSendOverlap = newMesh->getSendOverlap();
  assert(!newSendOverlap.isNull());
  Obj<mesh_type::recv_overlap_type> newRecvOverlap = newMesh->getRecvOverlap();
  assert(!newRecvOverlap.isNull());
  const Obj<mesh_type::send_overlap_type>& oldSendOverlap = oldMesh->getSendOverlap();
  assert(!oldSendOverlap.isNull());
  const Obj<mesh_type::recv_overlap_type>& oldRecvOverlap = oldMesh->getRecvOverlap();
  assert(!oldRecvOverlap.isNull());


  // Check edges in edgeToVertex for both endpoints sent to same process
  //   Put it in section with point being the lowest numbered vertex and value (other endpoint, new vertex)
  Obj<ALE::Section<point_type, EdgeType> > newVerticesSection = new ALE::Section<point_type, EdgeType>(oldMesh->comm());
  assert(!newVerticesSection.isNull());
  std::map<EdgeType, std::vector<int> > bndryEdgeToRank;
  
  const int localOffset = orderNewMesh.verticesNormal().min() - orderOldMesh.verticesNormal().min();

  for(std::map<EdgeType, point_type>::const_iterator e_iter = _edgeToVertex.begin(); e_iter != _edgeToVertex.end(); ++e_iter) {
    const point_type left  = e_iter->first.first;
    const point_type right = e_iter->first.second;
    
    if (oldSendOverlap->capContains(left) && oldSendOverlap->capContains(right)) {
      const Obj<mesh_type::send_overlap_type::traits::supportSequence>& leftRanksSeq = oldSendOverlap->support(left);
      assert(!leftRanksSeq.isNull());
      std::list<int> leftRanks(leftRanksSeq->begin(), leftRanksSeq->end());
      const Obj<mesh_type::send_overlap_type::traits::supportSequence>& rightRanks   = oldSendOverlap->support(right);
      assert(!rightRanks.isNull());
      std::list<int> ranks;
      std::set_intersection(leftRanks.begin(), leftRanks.end(), rightRanks->begin(), rightRanks->end(),
			    std::insert_iterator<std::list<int> >(ranks, ranks.begin()));
      
      if(ranks.size()) {
	newVerticesSection->addFiberDimension(std::min(e_iter->first.first, e_iter->first.second)+localOffset, 1);
	for(std::list<int>::const_iterator r_iter = ranks.begin(); r_iter != ranks.end(); ++r_iter) {
	  bndryEdgeToRank[e_iter->first].push_back(*r_iter);
	} // for
      } // if
    } // if
  } // for
  newVerticesSection->allocatePoint();
  const ALE::Section<point_type, EdgeType>::chart_type& chart = newVerticesSection->getChart();
  
  for(ALE::Section<point_type, EdgeType>::chart_type::const_iterator c_iter = chart.begin(); c_iter != chart.end(); ++c_iter) {
    typedef ALE::Section<point_type, EdgeType>::value_type value_type;
    const point_type p      = *c_iter;
    const int        dim    = newVerticesSection->getFiberDimension(p);
    int              v      = 0;
    value_type* values = (dim > 0) ? new value_type[dim] : 0;
    
    for(std::map<EdgeType, std::vector<int> >::const_iterator e_iter = bndryEdgeToRank.begin(); e_iter != bndryEdgeToRank.end() && v < dim; ++e_iter) {
      if (std::min(e_iter->first.first, e_iter->first.second)+localOffset == p) {
	values[v++] = EdgeType(std::max(e_iter->first.first, e_iter->first.second)+localOffset, _edgeToVertex[e_iter->first]);
      } // if
    } // for
    newVerticesSection->updatePoint(p, values);
    delete [] values;
  } // for
  // Copy across overlap
  typedef ALE::Pair<int, point_type> overlap_point_type;
  Obj<ALE::Section<overlap_point_type, EdgeType> > overlapVertices = new ALE::Section<overlap_point_type, EdgeType>(oldMesh->comm());
  
  ALE::Pullback::SimpleCopy::copy(newSendOverlap, newRecvOverlap, newVerticesSection, overlapVertices);
  // Merge by translating edge to local points, finding edge in _edgeToVertex, and adding (local new vetex, remote new vertex) to overlap
  for(std::map<EdgeType, std::vector<int> >::const_iterator e_iter = bndryEdgeToRank.begin(); e_iter != bndryEdgeToRank.end(); ++e_iter) {
    const point_type localPoint = _edgeToVertex[e_iter->first];
    
    for(std::vector<int>::const_iterator r_iter = e_iter->second.begin(); r_iter != e_iter->second.end(); ++r_iter) {
      point_type remoteLeft = -1, remoteRight = -1;
      const int  rank       = *r_iter;
      
      const Obj<mesh_type::send_overlap_type::traits::supportSequence>& leftRanks = newSendOverlap->support(e_iter->first.first+localOffset);
      for(mesh_type::send_overlap_type::traits::supportSequence::iterator lr_iter = leftRanks->begin(); lr_iter != leftRanks->end(); ++lr_iter) {
	if (rank == *lr_iter) {
	  remoteLeft = lr_iter.color();
	  break;
	} // if
      } // for
      const Obj<mesh_type::send_overlap_type::traits::supportSequence>& rightRanks = newSendOverlap->support(e_iter->first.second+localOffset);
      for(mesh_type::send_overlap_type::traits::supportSequence::iterator rr_iter = rightRanks->begin(); rr_iter != rightRanks->end(); ++rr_iter) {
	if (rank == *rr_iter) {
	  remoteRight = rr_iter.color();
	  break;
	} // if
      } // for
      const point_type remoteMin   = std::min(remoteLeft, remoteRight);
      const point_type remoteMax   = std::max(remoteLeft, remoteRight);
      const int        remoteSize  = overlapVertices->getFiberDimension(overlap_point_type(rank, remoteMin));
      const EdgeType *remoteVals  = overlapVertices->restrictPoint(overlap_point_type(rank, remoteMin));
      point_type       remotePoint = -1;
      
      for(int d = 0; d < remoteSize; ++d) {
	if (remoteVals[d].first == remoteMax) {
	  remotePoint = remoteVals[d].second;
	  break;
	} // if
      } // for
      newSendOverlap->addArrow(localPoint, rank, remotePoint);
      newRecvOverlap->addArrow(rank, localPoint, remotePoint);
    } // for
  } // for

  oldSendOverlap->view("OLD SEND OVERLAP");
  oldRecvOverlap->view("OLD RECV OVERLAP");
  newSendOverlap->view("NEW SEND OVERLAP");
  newRecvOverlap->view("NEW RECV OVERLAP");
} // overlapAddNewVertices

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
  static EdgeType splitEdges[6];
  
  assert(coneSize == 4);
  splitEdges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
  splitEdges[1] = EdgeType(std::min(cone[1], cone[2]), std::max(cone[1], cone[2]));
  splitEdges[2] = EdgeType(std::min(cone[2], cone[0]), std::max(cone[2], cone[0]));
  splitEdges[3] = EdgeType(std::min(cone[0], cone[3]), std::max(cone[0], cone[3]));
  splitEdges[4] = EdgeType(std::min(cone[1], cone[3]), std::max(cone[1], cone[3]));
  splitEdges[5] = EdgeType(std::min(cone[2], cone[3]), std::max(cone[2], cone[3]));
  *numEdges = 6;
  *edges    = splitEdges;
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
    static EdgeType splitEdges[9];

    assert(coneSize == 9);
    splitEdges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
    splitEdges[1] = EdgeType(std::min(cone[1], cone[2]), std::max(cone[1], cone[2]));
    splitEdges[2] = EdgeType(std::min(cone[2], cone[0]), std::max(cone[2], cone[0]));
    splitEdges[3] = EdgeType(std::min(cone[3], cone[4]), std::max(cone[3], cone[4]));
    splitEdges[4] = EdgeType(std::min(cone[4], cone[5]), std::max(cone[4], cone[5]));
    splitEdges[5] = EdgeType(std::min(cone[5], cone[3]), std::max(cone[5], cone[3]));
    splitEdges[6] = EdgeType(std::min(cone[6], cone[7]), std::max(cone[6], cone[7]));
    splitEdges[7] = EdgeType(std::min(cone[7], cone[8]), std::max(cone[7], cone[8]));
    splitEdges[8] = EdgeType(std::min(cone[8], cone[6]), std::max(cone[8], cone[6]));
    *numEdges = 9;
    *edges = splitEdges;
  } else {
    // Omit edges with censored (Lagrange multipler) vertices.
    static EdgeType splitEdges[6];

    assert(coneSize == 9);
    splitEdges[0] = EdgeType(std::min(cone[0], cone[1]), std::max(cone[0], cone[1]));
    splitEdges[1] = EdgeType(std::min(cone[1], cone[2]), std::max(cone[1], cone[2]));
    splitEdges[2] = EdgeType(std::min(cone[2], cone[0]), std::max(cone[2], cone[0]));
    splitEdges[3] = EdgeType(std::min(cone[3], cone[4]), std::max(cone[3], cone[4]));
    splitEdges[4] = EdgeType(std::min(cone[4], cone[5]), std::max(cone[4], cone[5]));
    splitEdges[5] = EdgeType(std::min(cone[5], cone[3]), std::max(cone[5], cone[3]));
    *numEdges = 6;
    *edges = splitEdges;
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
  const EdgeType  *edges;
  _edges_TETRAHEDRON(&edges, &numEdges, cone, coneSize);
  assert(numEdgesTet4 == numEdges);

  static point_type newCells[numNewCells*coneSizeTet4];
  point_type newVertices[numNewVertices];
  for(int iEdge=0, iNewVertex=0; iEdge < numEdgesTet4; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      throw ALE::Exception("Missing edge in refined mesh");
    } // if
    newVertices[iNewVertex++] = _edgeToVertex[edges[iEdge]];
  } // for

  // new cell 0
  newCells[0*4+0] = cone[0]+coneVertexOffset;
  newCells[0*4+1] = newVertices[3];
  newCells[0*4+2] = newVertices[0];
  newCells[0*4+3] = newVertices[2];

  // new cell 1
  newCells[1*4+0] = newVertices[0];
  newCells[1*4+1] = newVertices[1];
  newCells[1*4+2] = newVertices[2];
  newCells[1*4+3] = newVertices[3];

  // new cell 2
  newCells[2*4+0] = newVertices[0];
  newCells[2*4+1] = newVertices[3];
  newCells[2*4+2] = newVertices[4];
  newCells[2*4+3] = newVertices[1];

  // new cell 3
  newCells[3*4+0] = cone[1]+coneVertexOffset;
  newCells[3*4+1] = newVertices[4];
  newCells[3*4+2] = newVertices[1];
  newCells[3*4+3] = newVertices[0];

  // new cell 4
  newCells[4*4+0] = newVertices[2];
  newCells[4*4+1] = newVertices[5];
  newCells[4*4+2] = newVertices[3];
  newCells[4*4+3] = newVertices[1];

  // new cell 5
  newCells[5*4+0] = cone[2]+coneVertexOffset;
  newCells[5*4+1] = newVertices[5];
  newCells[5*4+2] = newVertices[2];
  newCells[5*4+3] = newVertices[1];

  // new cell 6
  newCells[6*4+0] = newVertices[1];
  newCells[6*4+1] = newVertices[4];
  newCells[6*4+2] = newVertices[5];
  newCells[6*4+3] = newVertices[3];

  // new cell 7
  newCells[7*4+0] = cone[3]+coneVertexOffset;
  newCells[7*4+1] = newVertices[3];
  newCells[7*4+2] = newVertices[5];
  newCells[7*4+3] = newVertices[4];

  *numCells = numNewCells;
  *cells = newCells;
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

  static point_type newCells[numNewCells*coneSizeTriPrism9];
  point_type newVertices[numNewVertices];
  for(int iEdge=0, iNewVertex=0; iEdge < numEdgesTriPrism9; ++iEdge) {
    if (_edgeToVertex.find(edges[iEdge]) == _edgeToVertex.end()) {
      throw ALE::Exception("Missing edge in refined mesh");
    } // if
    newVertices[iNewVertex++] = _edgeToVertex[edges[iEdge]];
  } // for

  newCells[0*9+0] = cone[0]+coneVertexOffsetNormal; // New cell 0
  newCells[0*9+1] = newVertices[0];
  newCells[0*9+2] = newVertices[2];
  newCells[0*9+3] = cone[3]+coneVertexOffsetNormal;
  newCells[0*9+4] = newVertices[3];
  newCells[0*9+5] = newVertices[5];
  newCells[0*9+6] = cone[6]+coneVertexOffsetCensored;
  newCells[0*9+7] = newVertices[6];
  newCells[0*9+8] = newVertices[8];
  
  newCells[1*9+0] = newVertices[0]; // New cell 1
  newCells[1*9+1] = newVertices[1];
  newCells[1*9+2] = newVertices[2];
  newCells[1*9+3] = newVertices[3];
  newCells[1*9+4] = newVertices[4];
  newCells[1*9+5] = newVertices[5];
  newCells[1*9+6] = newVertices[6];
  newCells[1*9+7] = newVertices[7];
  newCells[1*9+8] = newVertices[8];
  
  newCells[2*9+0] = cone[1]+coneVertexOffsetNormal; // New cell 2
  newCells[2*9+1] = newVertices[1];
  newCells[2*9+2] = newVertices[0];
  newCells[2*9+3] = cone[4]+coneVertexOffsetNormal;
  newCells[2*9+4] = newVertices[4];
  newCells[2*9+5] = newVertices[3];
  newCells[2*9+6] = cone[7]+coneVertexOffsetCensored;
  newCells[2*9+7] = newVertices[7];
  newCells[2*9+8] = newVertices[6];
  
  newCells[3*9+0] = cone[2]+coneVertexOffsetNormal; // New cell 3
  newCells[3*9+1] = newVertices[2];
  newCells[3*9+2] = newVertices[1];
  newCells[3*9+3] = cone[5]+coneVertexOffsetNormal;
  newCells[3*9+4] = newVertices[5];
  newCells[3*9+5] = newVertices[4];
  newCells[3*9+6] = cone[8]+coneVertexOffsetCensored;
  newCells[3*9+7] = newVertices[8];
  newCells[3*9+8] = newVertices[7];

  *numCells = numNewCells;
  *cells    = newCells;
} // _newCells_TRIANGLE_COHESIVE_LAGRANGE


// End of file 
