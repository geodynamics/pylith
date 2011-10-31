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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "RefineVol8Face4Edges2.hh" // implementation of class methods

#include "MeshOrder.hh" // USES MeshOrder

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
ALE::RefineVol8Face4Edges2::RefineVol8Face4Edges2(const mesh_type& mesh) :
  _mesh(mesh)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
ALE::RefineVol8Face4Edges2::~RefineVol8Face4Edges2(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set coordinates of new vertices.
void
ALE::RefineVol8Face4Edges2::setCoordsNewVertices(const ALE::Obj<mesh_type::real_section_type>& newCoordsSection,
						 const ALE::Obj<mesh_type::real_section_type>& oldCoordsSection)
{ // setCoordsNewVertices
  assert(!newCoordsSection.isNull());
  assert(!oldCoordsSection.isNull());

  PylithScalar coordinatesVertex[3];

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

  const face_map_type::const_iterator facesEnd = _faceToVertex.end();
  for (face_map_type::const_iterator f_iter = _faceToVertex.begin(); f_iter != facesEnd; ++f_iter) {
    const point_type newVertex = f_iter->second;

    assert(spaceDim == newCoordsSection->getFiberDimension(newVertex));
    for (int iDim=0; iDim < spaceDim; ++iDim)
      coordinatesVertex[iDim] = 0.0;
    for (int iVertex=0; iVertex < 4; ++iVertex) {
      const point_type faceVertex = f_iter->first.points[iVertex];
      assert(spaceDim == oldCoordsSection->getFiberDimension(faceVertex));

      const mesh_type::real_section_type::value_type* coords = oldCoordsSection->restrictPoint(faceVertex);
      for (int iDim=0; iDim < spaceDim; ++iDim)
	coordinatesVertex[iDim] += 0.25*coords[iDim];
    } // for
    newCoordsSection->updatePoint(newVertex, coordinatesVertex);
  } // for

  const volume_map_type::const_iterator volumesEnd = _volumeToVertex.end();
  for (volume_map_type::const_iterator v_iter = _volumeToVertex.begin(); v_iter != volumesEnd; ++v_iter) {
    const point_type newVertex = v_iter->second;

    assert(spaceDim == newCoordsSection->getFiberDimension(newVertex));
    for (int iDim=0; iDim < spaceDim; ++iDim)
      coordinatesVertex[iDim] = 0.0;
    for (int iVertex=0; iVertex < 8; ++iVertex) {
      const point_type volumeVertex = v_iter->first.points[iVertex];
      assert(spaceDim == oldCoordsSection->getFiberDimension(volumeVertex));

      const mesh_type::real_section_type::value_type* coords = oldCoordsSection->restrictPoint(volumeVertex);
      for (int iDim=0; iDim < spaceDim; ++iDim)
	coordinatesVertex[iDim] += 0.125*coords[iDim];
    } // for
    newCoordsSection->updatePoint(newVertex, coordinatesVertex);
  } // for
} // setCoordsNewVertices

// ----------------------------------------------------------------------
// Add space for new vertices in group.
void
ALE::RefineVol8Face4Edges2::groupAddNewVertices(const ALE::Obj<mesh_type::int_section_type>& newGroup,
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

  const face_map_type::const_iterator facesEnd = _faceToVertex.end();
  for (face_map_type::const_iterator f_iter = _faceToVertex.begin(); f_iter != facesEnd; ++f_iter) {
    const point_type newVertex = f_iter->second;

    bool hasFace = true;
    for (int iVertex=0; iVertex < 4; ++iVertex) {
      const point_type faceVertex = f_iter->first.points[iVertex];
      if (!oldGroup->getFiberDimension(faceVertex)) {
	hasFace = false;
	break;
      } // if
    } // for
    if (hasFace) {
      newGroup->setFiberDimension(newVertex, 1);
    } // if
  } // for

  const volume_map_type::const_iterator volumesEnd = _volumeToVertex.end();
  for (volume_map_type::const_iterator v_iter = _volumeToVertex.begin(); v_iter != volumesEnd; ++v_iter) {
    const point_type newVertex = v_iter->second;

    bool hasVolume = true;
    for (int iVertex=0; iVertex < 8; ++iVertex) {
      const point_type volumeVertex = v_iter->first.points[iVertex];
      if (!oldGroup->getFiberDimension(volumeVertex)) {
	hasVolume = false;
	break;
      } // if
    } // for
    if (hasVolume) {
      newGroup->setFiberDimension(newVertex, 1);
    } // if
  } // for
} // groupAddNewVertices

// ----------------------------------------------------------------------
// Set new vertices in group.
void
ALE::RefineVol8Face4Edges2::groupSetNewVertices(const ALE::Obj<mesh_type::int_section_type>& newGroup,
						const ALE::Obj<mesh_type::int_section_type>& oldGroup)
{ // groupSetNewVertices
  assert(!newGroup.isNull());
  assert(!oldGroup.isNull());

  const edge_map_type::const_iterator edgesEnd = _edgeToVertex.end();
  for (edge_map_type::const_iterator e_iter = _edgeToVertex.begin(); e_iter != edgesEnd; ++e_iter) {
    const point_type newVertex = e_iter->second;
    const point_type edgeVertex = e_iter->first.first;

    if (1 == newGroup->getFiberDimension(newVertex)) {
      newGroup->updatePoint(newVertex, oldGroup->restrictPoint(edgeVertex));
      //std::cout << "Adding new vertex: " << newVertex << " based on edge " << e_iter->first << std::endl;
    } // if
  } // for

  const face_map_type::const_iterator facesEnd = _faceToVertex.end();
  for (face_map_type::const_iterator f_iter = _faceToVertex.begin(); f_iter != facesEnd; ++f_iter) {
    const point_type newVertex = f_iter->second;
    const point_type faceVertex = f_iter->first.points[0];

    if (1 == newGroup->getFiberDimension(newVertex)) {
      newGroup->updatePoint(newVertex, oldGroup->restrictPoint(faceVertex));
      //std::cout << "Adding new vertex: " << newVertex << " based on face " << f_iter->first << std::endl;
    } // if
  } // for

  const volume_map_type::const_iterator volumesEnd = _volumeToVertex.end();
  for (volume_map_type::const_iterator v_iter = _volumeToVertex.begin(); v_iter != volumesEnd; ++v_iter) {
    const point_type newVertex = v_iter->second;
    const point_type volumeVertex = v_iter->first.points[0];

    if (1 == newGroup->getFiberDimension(newVertex)) {
      newGroup->updatePoint(newVertex, oldGroup->restrictPoint(volumeVertex));
      //std::cout << "Adding new vertex: " << newVertex << " based on volume " << v_iter->first << std::endl;
    } // if
  } // for
} // groupSetNewVertices

// ----------------------------------------------------------------------
// Add new vertices to label.
void
ALE::RefineVol8Face4Edges2::labelAddNewVertices(const ALE::Obj<mesh_type>& newMesh,
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

  const face_map_type::const_iterator facesEnd = _faceToVertex.end();
  for (face_map_type::const_iterator f_iter = _faceToVertex.begin(); f_iter != facesEnd; ++f_iter) {
    const point_type newVertex = f_iter->second;
    const point_type faceVertex = f_iter->first.points[0];
    const int value = oldMesh->getValue(oldLabel, faceVertex, defaultValue);

    if (value != defaultValue) {
      bool hasFace = true;
      for (int iVertex=1; iVertex < 4; ++iVertex) {
	const point_type faceVertex2 = f_iter->first.points[iVertex];
	const int value2 = oldMesh->getValue(oldLabel, faceVertex2, defaultValue);
	if (value2 != value) {
	  hasFace = false;
	  break;
	} // if
      } // for
      if (hasFace) {
	newMesh->setValue(newLabel, newVertex, value);
      } // if
    } // if
  } // for

  const volume_map_type::const_iterator volumesEnd = _volumeToVertex.end();
  for (volume_map_type::const_iterator v_iter = _volumeToVertex.begin(); v_iter != volumesEnd; ++v_iter) {
    const point_type newVertex = v_iter->second;
    const point_type volumeVertex = v_iter->first.points[0];
    const int value = oldMesh->getValue(oldLabel, volumeVertex, defaultValue);

    if (value != defaultValue) {
      bool hasVolume = true;
      for (int iVertex=1; iVertex < 8; ++iVertex) {
	const point_type volumeVertex2 = v_iter->first.points[iVertex];
	const int value2 = oldMesh->getValue(oldLabel, volumeVertex2, defaultValue);
	if (value2 != value) {
	  hasVolume = false;
	  break;
	} // if
      } // for
      if (hasVolume) {
	newMesh->setValue(newLabel, newVertex, value);
      } // if
    } // if
  } // for
} // labelAddNewVertices

// ----------------------------------------------------------------------
// Calculate new overlap.
void
ALE::RefineVol8Face4Edges2::overlapAddNewVertices(const Obj<mesh_type>& newMesh,
					    const MeshOrder& orderNewMesh,
					    const Obj<mesh_type>& oldMesh,
					    const MeshOrder& orderOldMesh)
{ // overlapAddNewVertices
  assert(!newMesh.isNull());
  assert(!oldMesh.isNull());

  // :TODO: Add face vertices

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
  //     Notice that points are converted to the new numbering with refined cells
  Obj<ALE::Section<point_type, EdgeType> > newVerticesSection = new ALE::Section<point_type, EdgeType>(oldMesh->comm());
  assert(!newVerticesSection.isNull());
  std::map<EdgeType, std::vector<int> > bndryEdgeToRank; // Maps an edge to a set of ranks who share it
  
  const int localNormalOffset   = orderNewMesh.verticesNormal().min()   - orderOldMesh.verticesNormal().min();
  const int localCensoredOffset = orderNewMesh.verticesCensored().min() - orderOldMesh.verticesCensored().min();

  for(std::map<EdgeType, point_type>::const_iterator e_iter = _edgeToVertex.begin(); e_iter != _edgeToVertex.end(); ++e_iter) {
    const point_type left  = e_iter->first.first;
    const point_type right = e_iter->first.second;
    
    if (oldSendOverlap->capContains(left) && oldSendOverlap->capContains(right)) {
      const Obj<mesh_type::send_overlap_type::supportSequence>& leftRanksSeq  = oldSendOverlap->support(left);
      assert(!leftRanksSeq.isNull());
      std::set<int> leftRanks(leftRanksSeq->begin(), leftRanksSeq->end());
      const Obj<mesh_type::send_overlap_type::supportSequence>& rightRanksSeq = oldSendOverlap->support(right);
      assert(!rightRanksSeq.isNull());
      std::set<int> rightRanks(rightRanksSeq->begin(), rightRanksSeq->end());
      std::set<int> ranks;
      std::set_intersection(leftRanks.begin(), leftRanks.end(), rightRanks.begin(), rightRanks.end(),
			    std::insert_iterator<std::set<int> >(ranks, ranks.begin()));
      
      if(ranks.size()) {
        const int localOffset = orderOldMesh.verticesNormal().hasPoint(e_iter->first.first) ? localNormalOffset : localCensoredOffset;

        assert(orderOldMesh.verticesNormal().hasPoint(e_iter->first.first) == orderOldMesh.verticesNormal().hasPoint(e_iter->first.second));
        newVerticesSection->addFiberDimension(std::min(e_iter->first.first, e_iter->first.second)+localOffset, 1);
        for(std::set<int>::const_iterator r_iter = ranks.begin(); r_iter != ranks.end(); ++r_iter) {
          bndryEdgeToRank[e_iter->first].push_back(*r_iter);
        } // for
      } // if
    } // if
  } // for
  newVerticesSection->allocatePoint();
  const ALE::Section<point_type, EdgeType>::chart_type& edgeChart = newVerticesSection->getChart();
  
  for(ALE::Section<point_type, EdgeType>::chart_type::const_iterator c_iter = edgeChart.begin(); c_iter != edgeChart.end(); ++c_iter) {
    typedef ALE::Section<point_type, EdgeType>::value_type value_type;
    const point_type p      = *c_iter;
    const int        dim    = newVerticesSection->getFiberDimension(p);
    int              v      = 0;
    value_type* values = (dim > 0) ? new value_type[dim] : 0;
    
    for(std::map<EdgeType, std::vector<int> >::const_iterator e_iter = bndryEdgeToRank.begin(); e_iter != bndryEdgeToRank.end() && v < dim; ++e_iter) {
      const int localOffset = orderOldMesh.verticesNormal().hasPoint(e_iter->first.first) ? localNormalOffset : localCensoredOffset;

      assert(orderOldMesh.verticesNormal().hasPoint(e_iter->first.first) == orderOldMesh.verticesNormal().hasPoint(e_iter->first.second));
      if (std::min(e_iter->first.first, e_iter->first.second)+localOffset == p) {
        values[v++] = EdgeType(std::max(e_iter->first.first, e_iter->first.second)+localOffset, _edgeToVertex[e_iter->first]);
      } // if
    } // for
    newVerticesSection->updatePoint(p, values);
    delete [] values;
  } // for

  // Check faces in faceToVertex for all corners sent to same process
  //   Put it in section with point being the lowest numbered vertex and value (other endpoints, new vertex)
  //     Notice that points are converted to the new numbering with refined cells
  Obj<ALE::Section<point_type, FaceType> > newFaceVerticesSection = new ALE::Section<point_type, FaceType>(oldMesh->comm());
  assert(!newFaceVerticesSection.isNull());
  std::map<FaceType, std::vector<int>, FaceCmp<point_type>  > bndryFaceToRank;

  for(std::map<FaceType, point_type>::const_iterator f_iter = _faceToVertex.begin(); f_iter != _faceToVertex.end(); ++f_iter) {
    bool isRemote = true;
    
    for(int i = 0; i < 4; ++i) {
      if (!oldSendOverlap->capContains(f_iter->first.points[i])) {
        isRemote = false;
        break;
      }
    }
    if (isRemote) {
      const point_type first     = f_iter->first.points[0];
      point_type       minVertex = first;

      const Obj<mesh_type::send_overlap_type::supportSequence>& firstRanksSeq  = oldSendOverlap->support(first);
      assert(!firstRanksSeq.isNull());
      std::set<int> firstRanks(firstRanksSeq->begin(), firstRanksSeq->end());
      std::set<int> ranks, curRanks;

      curRanks = firstRanks;
      for(int i = 1; i < 4; ++i) {
        const point_type nextVertex = f_iter->first.points[i];

        const Obj<mesh_type::send_overlap_type::supportSequence>& nextRanksSeq = oldSendOverlap->support(nextVertex);
        assert(!nextRanksSeq.isNull());
        std::set<int> nextRanks(nextRanksSeq->begin(), nextRanksSeq->end());

        ranks.clear();
        std::set_intersection(curRanks.begin(), curRanks.end(), nextRanks.begin(), nextRanks.end(),
                              std::insert_iterator<std::set<int> >(ranks, ranks.begin()));
        curRanks = ranks;
        minVertex = std::min(nextVertex, minVertex);
      }
      if (ranks.size()) {
        const int localOffset = orderOldMesh.verticesNormal().hasPoint(minVertex) ? localNormalOffset : localCensoredOffset;

        assert(orderNewMesh.verticesNormal().hasPoint(f_iter->second));
        newFaceVerticesSection->addFiberDimension(minVertex+localOffset, 1);
        for(std::set<int>::const_iterator r_iter = ranks.begin(); r_iter != ranks.end(); ++r_iter) {
          bndryFaceToRank[f_iter->first].push_back(*r_iter);
        }
      }
    }
  }
  newFaceVerticesSection->allocatePoint();
  const ALE::Section<point_type, FaceType>::chart_type& faceChart = newFaceVerticesSection->getChart();
  
  for(ALE::Section<point_type, FaceType>::chart_type::const_iterator c_iter = faceChart.begin(); c_iter != faceChart.end(); ++c_iter) {
    typedef ALE::Section<point_type, FaceType>::value_type value_type;
    const point_type p      = *c_iter;
    const int        dim    = newFaceVerticesSection->getFiberDimension(p);
    int              v      = 0;
    value_type* values = (dim > 0) ? new value_type[dim] : 0;
    
    for(std::map<FaceType, std::vector<int> >::const_iterator f_iter = bndryFaceToRank.begin(); f_iter != bndryFaceToRank.end() && v < dim; ++f_iter) {
      const point_type first     = f_iter->first.points[0];
      point_type       minVertex = first;
      int              localOffset;

      for(int i = 1; i < 4; ++i) {
        const point_type nextVertex = f_iter->first.points[i];
        minVertex = std::min(nextVertex, minVertex);
      }
      localOffset = orderOldMesh.verticesNormal().hasPoint(minVertex) ? localNormalOffset : localCensoredOffset;

      if (minVertex+localOffset == p) {
        FaceType face;
        int      k = 0;

        for(int i = 0; i < 4; ++i) {
          if (f_iter->first.points[i] != minVertex) {
            localOffset = orderOldMesh.verticesNormal().hasPoint(f_iter->first.points[i]) ? localNormalOffset : localCensoredOffset;

            face.points[k++] = f_iter->first.points[i]+localOffset;
          }
        }
        assert(k == 3);
        std::sort(&face.points[0], &face.points[3]);
        face.points[3] = _faceToVertex[f_iter->first];
        values[v++] = face;
      } // if
    } // for
    assert(v == dim);
    newFaceVerticesSection->updatePoint(p, values);
    delete [] values;
  } // for

  // Copy across overlap
  typedef ALE::Pair<int, point_type> overlap_point_type;
  // This maps (remote rank, remote min point) --> (remote max point, new remote vertex)
  Obj<ALE::Section<overlap_point_type, EdgeType> > overlapVertices     = new ALE::Section<overlap_point_type, EdgeType>(oldMesh->comm());
  // This maps (remote rank, remote min point) --> (other remote points..., new remote vertex)
  Obj<ALE::Section<overlap_point_type, FaceType> > overlapFaceVertices = new ALE::Section<overlap_point_type, FaceType>(oldMesh->comm());
  
  ALE::Pullback::SimpleCopy::copy(newSendOverlap, newRecvOverlap, newVerticesSection,     overlapVertices);
  ALE::Pullback::SimpleCopy::copy(newSendOverlap, newRecvOverlap, newFaceVerticesSection, overlapFaceVertices);
  // Merge by translating edge to local points, finding edge in _edgeToVertex, and adding (local new vetex, remote new vertex) to overlap
  //   Loop over all shared edges
  for(std::map<EdgeType, std::vector<int> >::const_iterator e_iter = bndryEdgeToRank.begin(); e_iter != bndryEdgeToRank.end(); ++e_iter) {
    // Point added on this edge by refinement
    const point_type localPoint = _edgeToVertex[e_iter->first];

    // Loop over all ranks which share this edge
    for(std::vector<int>::const_iterator r_iter = e_iter->second.begin(); r_iter != e_iter->second.end(); ++r_iter) {
      point_type remoteLeft   = -1, remoteRight = -1;
      const int  rank         = *r_iter;
      // Offset for left edge point in the new mesh (check whether edge is between Lagrange vertices)
      const int  localOffsetL = orderOldMesh.verticesNormal().hasPoint(e_iter->first.first) ? localNormalOffset : localCensoredOffset;

      // Find the rank which owns the left edge point, and return the remotePoint in the new mesh
      const Obj<mesh_type::send_overlap_type::supportSequence>& leftRanks = newSendOverlap->support(e_iter->first.first+localOffsetL);
      for(mesh_type::send_overlap_type::supportSequence::iterator lr_iter = leftRanks->begin(); lr_iter != leftRanks->end(); ++lr_iter) {
        if (rank == *lr_iter) {
          remoteLeft = lr_iter.color();
          break;
        } // if
      } // for
      // Offset for right edge point in the new mesh (check whether edge is between Lagrange vertices)
      const int  localOffsetR = orderOldMesh.verticesNormal().hasPoint(e_iter->first.second) ? localNormalOffset : localCensoredOffset;

      // Find the rank which owns the right edge point, and return the remotePoint in the new mesh
      const Obj<mesh_type::send_overlap_type::supportSequence>& rightRanks = newSendOverlap->support(e_iter->first.second+localOffsetR);
      for(mesh_type::send_overlap_type::supportSequence::iterator rr_iter = rightRanks->begin(); rr_iter != rightRanks->end(); ++rr_iter) {
        if (rank == *rr_iter) {
          remoteRight = rr_iter.color();
          break;
        } // if
      } // for
      const point_type remoteMin   = std::min(remoteLeft, remoteRight);
      const point_type remoteMax   = std::max(remoteLeft, remoteRight);
      const int        remoteSize  = overlapVertices->getFiberDimension(overlap_point_type(rank, remoteMin));
      const EdgeType  *remoteVals  = overlapVertices->restrictPoint(overlap_point_type(rank, remoteMin));
      point_type       remotePoint = -1;
      
      // Loop over all shared edges from rank with endpoint remoteMin
      for(int d = 0; d < remoteSize; ++d) {
        if (remoteVals[d].first == remoteMax) {
          remotePoint = remoteVals[d].second;
          break;
        } // if
      } // for
#if 0
      // Debugging for fault edge problem
      if (remotePoint < 0) {
        std::cout << "["<<oldMesh->commRank()<<"]Failed to find remote partner for edge " << e_iter->first << " from rank " << rank << std::endl;
        std::cout << "["<<oldMesh->commRank()<<"]   remote edge ("<<remoteLeft<<","<<remoteRight<<") had remoteSize " << remoteSize << std::endl;
        for(int d = 0; d < remoteSize; ++d) {
          std::cout << "["<<oldMesh->commRank()<<"]     remote val " << remoteVals[d] << std::endl;
        }
        const Obj<mesh_type::send_overlap_type::supportSequence>& leftRanks2 = oldSendOverlap->support(e_iter->first.first);
        for(mesh_type::send_overlap_type::supportSequence::iterator lr_iter = leftRanks2->begin(); lr_iter != leftRanks2->end(); ++lr_iter) {
          if (rank == *lr_iter) {
            std::cout << "["<<oldMesh->commRank()<<"]     left match:  old vertex " << lr_iter.color() << std::endl;
            break;
          }
        }
        const Obj<mesh_type::send_overlap_type::supportSequence>& rightRanks2 = oldSendOverlap->support(e_iter->first.second);
        for(mesh_type::send_overlap_type::supportSequence::iterator rr_iter = rightRanks2->begin(); rr_iter != rightRanks2->end(); ++rr_iter) {
          if (rank == *rr_iter) {
            std::cout << "["<<oldMesh->commRank()<<"]     right match: old vertex " << rr_iter.color() << std::endl;
            break;
          }
        }
      }
      assert(remotePoint >= 0);
#endif
      if (remotePoint >= 0) {
        newSendOverlap->addArrow(localPoint, rank, remotePoint);
        newRecvOverlap->addArrow(rank, localPoint, remotePoint);
      }
    } // for
  } // for
  // Merge by translating face to local points, finding face in _faceToVertex, and adding (local new vetex, remote new vertex) to overlap
  for(std::map<FaceType, std::vector<int> >::const_iterator f_iter = bndryFaceToRank.begin(); f_iter != bndryFaceToRank.end(); ++f_iter) {
    const point_type localPoint = _faceToVertex[f_iter->first];
    
    for(std::vector<int>::const_iterator r_iter = f_iter->second.begin(); r_iter != f_iter->second.end(); ++r_iter) {
      FaceType  remoteVertices(-1); // These are the remote vertices on process 'rank' for this face
      const int rank = *r_iter;

      for(int i = 0; i < 4; ++i) {
        const int localOffset = orderOldMesh.verticesNormal().hasPoint(f_iter->first.points[i]) ? localNormalOffset : localCensoredOffset;
        const Obj<mesh_type::send_overlap_type::supportSequence>& faceRanks = newSendOverlap->support(f_iter->first.points[i]+localOffset);
        for(mesh_type::send_overlap_type::supportSequence::iterator fr_iter = faceRanks->begin(); fr_iter != faceRanks->end(); ++fr_iter) {
          if (rank == *fr_iter) {
            remoteVertices.points[i] = fr_iter.color();
            break;
          } // if
        } // for
        assert(remoteVertices.points[i] >= 0);
      }
      const point_type remoteMin   = std::min(std::min(std::min(remoteVertices.points[0], remoteVertices.points[1]), remoteVertices.points[2]), remoteVertices.points[3]);
      const int        remoteSize  = overlapFaceVertices->getFiberDimension(overlap_point_type(rank, remoteMin));
      const FaceType  *remoteVals  = overlapFaceVertices->restrictPoint(overlap_point_type(rank, remoteMin));
      point_type       remotePoint = -1;
      int              k           = 0;
      FaceType         remoteMax;

      for(int i = 0; i < 4; ++i) {
        if (remoteVertices.points[i] == remoteMin) continue;
        remoteMax.points[k++] = remoteVertices.points[i];
      }
      assert(k == 3);
      std::sort(&remoteMax.points[0], &remoteMax.points[3]);
      for(int d = 0; d < remoteSize; ++d) {
        int i = 0;

        for(i = 0; i < 3; ++i) {
          if (remoteVals[d].points[i] != remoteMax.points[i]) break;
        }
        if (i == 3) {
          remotePoint = remoteVals[d].points[3];
          break;
        } // if
      } // for
      assert(localPoint >= orderNewMesh.verticesNormal().min() && localPoint < orderNewMesh.verticesNormal().max());
#if 0
      // Debugging for fault edge problem
      assert(remotePoint >= 0);
#endif
      if (remotePoint >= 0) {
        newSendOverlap->addArrow(localPoint, rank, remotePoint);
        newRecvOverlap->addArrow(rank, localPoint, remotePoint);
      }
    } // for
  } // for

#if 0 // debuggin
  oldSendOverlap->view("OLD SEND OVERLAP");
  oldRecvOverlap->view("OLD RECV OVERLAP");
  newSendOverlap->view("NEW SEND OVERLAP");
  newRecvOverlap->view("NEW RECV OVERLAP");
#endif
} // overlapAddNewVertces


// End of file 
