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

#include "MeshRefiner.hh" // implementation of class methods

#include "MeshOrder.hh" // USES MeshOrder

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
template<typename cellrefiner_type>
ALE::MeshRefiner<cellrefiner_type>::MeshRefiner(void) :
  _orderOldMesh(new MeshOrder),
  _orderNewMesh(new MeshOrder)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename cellrefiner_type>
ALE::MeshRefiner<cellrefiner_type>::~MeshRefiner(void)
{ // destructor
  delete _orderOldMesh; _orderOldMesh = PETSC_NULL;
  delete _orderNewMesh; _orderNewMesh = PETSC_NULL;
} // destructor

// ----------------------------------------------------------------------
// Refine mesh.
template<typename cellrefiner_type>
void
ALE::MeshRefiner<cellrefiner_type>::refine(const Obj<mesh_type>& newMesh, 
				      const Obj<mesh_type>& mesh, 
				      cellrefiner_type& refiner)
{ // refine
  assert(!mesh.isNull());
  if (mesh->hasLabel("censored depth")) {
    _refineCensored(newMesh, mesh, refiner);
  } else {
    _refine(newMesh, mesh, refiner);
  } // if/else
} // refine

// ----------------------------------------------------------------------
// Refine mesh.
template<typename cellrefiner_type>
void
ALE::MeshRefiner<cellrefiner_type>::_refine(const Obj<mesh_type>& newMesh, 
				       const Obj<mesh_type>& mesh, 
				       cellrefiner_type& refiner)
{ // _refine
  typedef Interval<point_type> interval_type;

  assert(_orderOldMesh);
  assert(_orderNewMesh);

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //logger.setDebug(1);
  logger.stagePush("Mesh");
  logger.stagePush("MeshCreation");

  // Calculate order in old mesh.
  _orderOldMesh->initialize(mesh);

  assert(!mesh.isNull());
  assert(!newMesh.isNull());
  
  // Get original mesh stuff.
  const Obj<mesh_type::label_sequence>& cells = mesh->heightStratum(0);
  assert(!cells.isNull());
  const mesh_type::label_sequence::iterator cellsEnd = cells->end();
  
  const Obj<mesh_type::label_sequence>& vertices = mesh->depthStratum(0);
  assert(!vertices.isNull());
  
  const Obj<mesh_type::sieve_type>& sieve = mesh->getSieve();
  assert(!sieve.isNull());
  ALE::ISieveVisitor::PointRetriever<mesh_type::sieve_type> cV(std::max(1, sieve->getMaxConeSize()));
  
  // Count number of cells.
  int newNumCells = 0;
  for (mesh_type::label_sequence::iterator c_iter = cells->begin(); c_iter != cellsEnd; ++c_iter) {
    newNumCells += refiner.numNewCells(*c_iter);
  } // for

  // Count number of vertices.
  const int oldNumVertices = vertices->size();
  int counterBegin = newNumCells + oldNumVertices;
  point_type curNewVertex = counterBegin;
  for(mesh_type::label_sequence::iterator c_iter = cells->begin(); c_iter != cellsEnd; ++c_iter) {
    cV.clear();
    sieve->cone(*c_iter, cV);
    refiner.splitCell(*c_iter, cV.getPoints(), cV.getSize(), &curNewVertex);
  } // for
  const int newNumVertices = oldNumVertices + curNewVertex - counterBegin;

  _orderNewMesh->cellsNormal(0, newNumCells);
  _orderNewMesh->verticesNormal(newNumCells, newNumCells+newNumVertices);
  _orderNewMesh->verticesCensored(newNumCells+newNumVertices, newNumCells+newNumVertices);
  _orderNewMesh->cellsCensored(newNumCells+newNumVertices, newNumCells+newNumVertices);
  
  // Allocate chart for new sieve.
  const Obj<mesh_type::sieve_type>& newSieve = newMesh->getSieve();
  assert(!newSieve.isNull());
  newSieve->setChart(mesh_type::sieve_type::chart_type(0, _orderNewMesh->cellsCensored().max()));

  // Create new sieve with correct sizes for refined cells
  point_type curNewCell = _orderNewMesh->cellsNormal().min();
  const interval_type::const_iterator oldCellsEnd = _orderOldMesh->cellsNormal().end();
  for (interval_type::const_iterator c_iter=_orderOldMesh->cellsNormal().begin(); c_iter != oldCellsEnd; ++c_iter) {
    // Set new cone and support sizes
    cV.clear();
    sieve->cone(*c_iter, cV);
    const point_type* cone = cV.getPoints();
    const int coneSize = cV.getSize();

    const point_type* newCells;
    int numNewCells = 0;
    refiner.getNewCells(&newCells, &numNewCells, *c_iter, cone, coneSize, *_orderOldMesh, *_orderNewMesh);

    for(int iCell=0; iCell < numNewCells; ++iCell, ++curNewCell) {
      newSieve->setConeSize(curNewCell, coneSize);
      for(int iVertex=0; iVertex < coneSize; ++iVertex) {
	newSieve->addSupportSize(newCells[iCell*coneSize+iVertex], 1);
      } // for
    } // for

  } // for
  newSieve->allocate();

  // Create refined cells in new sieve.
  curNewCell = _orderNewMesh->cellsNormal().min();
  for (interval_type::const_iterator c_iter=_orderOldMesh->cellsNormal().begin(); c_iter != oldCellsEnd; ++c_iter) {
    cV.clear();
    sieve->cone(*c_iter, cV);
    const point_type *cone = cV.getPoints();
    const int coneSize = cV.getSize();
    
    const point_type* newCells;
    int numNewCells = 0;
    refiner.getNewCells(&newCells, &numNewCells, *c_iter, cone, coneSize, *_orderOldMesh, *_orderNewMesh);
    
    for(int iCell=0; iCell < numNewCells; ++iCell, ++curNewCell) {
      newSieve->setCone(&newCells[iCell*coneSize], curNewCell);
    } // for
  } // for
  newSieve->symmetrize();


  logger.stagePop();
  //logger.setDebug(0);
  logger.stagePush("MeshCoordinates");

  // Set coordinates in refined mesh.
  const Obj<mesh_type::real_section_type>& coordinates = mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const Obj<mesh_type::real_section_type>& newCoordinates = newMesh->getRealSection("coordinates");
  assert(!newCoordinates.isNull());

  const mesh_type::label_sequence::const_iterator verticesEnd = vertices->end();
  assert(vertices->size() > 0);
  const int spaceDim = coordinates->getFiberDimension(*vertices->begin());
  assert(spaceDim > 0);
  newCoordinates->setChart(mesh_type::sieve_type::chart_type(_orderNewMesh->verticesNormal().min(), _orderNewMesh->verticesCensored().max()));

  interval_type::const_iterator newVerticesEnd = _orderNewMesh->verticesCensored().end();
  for (interval_type::const_iterator v_iter=_orderNewMesh->verticesNormal().begin(); v_iter != newVerticesEnd; ++v_iter) {
    newCoordinates->setFiberDimension(*v_iter, spaceDim);
  } // for
  newCoordinates->allocatePoint();
      
  const interval_type::const_iterator oldVerticesEnd = _orderOldMesh->verticesCensored().end();
  for (interval_type::const_iterator vOld_iter=_orderOldMesh->verticesNormal().begin(), vNew_iter=_orderNewMesh->verticesNormal().begin(); vOld_iter != oldVerticesEnd; ++vOld_iter, ++vNew_iter) {
    //std::cout << "Copy coordinates from old vertex " << *vOld_iter << " to new vertex " << *vNew_iter << std::endl;
    newCoordinates->updatePoint(*vNew_iter, coordinates->restrictPoint(*vOld_iter));
  } // for

  refiner.setCoordsNewVertices(newCoordinates, coordinates);

  logger.stagePop();
  logger.stagePush("MeshStratification");

  _stratify(newMesh);

  logger.stagePop();
  logger.stagePush("MeshOverlap");

  _calcNewOverlap(newMesh, mesh, refiner);

  logger.stagePop();
  logger.stagePush("MeshIntSections");

  _createIntSections(newMesh, mesh, refiner);

  logger.stagePop();
  logger.stagePush("MeshLabels");

  _createLabels(newMesh, mesh, refiner);

  logger.stagePop();
  logger.stagePop(); // Mesh
} // _refine
  
// ----------------------------------------------------------------------
// Refine mesh with a censored depth.
template<typename cellrefiner_type>
void
ALE::MeshRefiner<cellrefiner_type>::_refineCensored(const Obj<mesh_type>& newMesh, 
					       const Obj<mesh_type>& mesh, 
					       cellrefiner_type& refiner)
{ // _refineCensored
  typedef Interval<point_type> interval_type;

  assert(_orderOldMesh);
  assert(_orderNewMesh);

  // Calculate order in old mesh.
  _orderOldMesh->initialize(mesh);

  assert(!mesh.isNull());
  assert(!newMesh.isNull());
  
  // Get original mesh stuff.
  const Obj<mesh_type::label_sequence>& cells = mesh->heightStratum(0);
  assert(!cells.isNull());
  const mesh_type::label_sequence::iterator cellsEnd = cells->end();
  
  const Obj<mesh_type::label_sequence>& vertices = mesh->depthStratum(0);
  assert(!vertices.isNull());
  
  const Obj<mesh_type::sieve_type>& sieve = mesh->getSieve();
  assert(!sieve.isNull());
  ALE::ISieveVisitor::PointRetriever<mesh_type::sieve_type> cV(std::max(1, sieve->getMaxConeSize()));

  // Count number of cells in censored depth (normal cells).
  const Obj<mesh_type::label_sequence>& cellsNormal = mesh->getLabelStratum("censored depth", mesh->depth());
  assert(!cellsNormal.isNull());
  const mesh_type::label_sequence::iterator cellsNormalEnd = cellsNormal->end();
  int newNumCellsNormal = 0;
  for(mesh_type::label_sequence::iterator c_iter = cellsNormal->begin(); c_iter != cellsNormalEnd; ++c_iter)
    newNumCellsNormal += refiner.numNewCells(*c_iter);
	
  // Count number of remaining cells (other cells).
  const int numSkip = cellsNormal->size();  
  mesh_type::label_sequence::const_iterator c_iter = cells->begin();
  for (int i=0; i < numSkip; ++i)
    ++c_iter;
  const mesh_type::label_sequence::const_iterator cellsCensoredBegin = c_iter;
  int newNumCellsCensored = 0;
  for (; c_iter != cellsEnd; ++c_iter)
    newNumCellsCensored += refiner.numNewCells(*c_iter);
  
  // Count number of normal vertices.
  assert(!mesh->getFactory().isNull());
  Obj<mesh_type::numbering_type> vNumbering = mesh->getFactory()->getNumbering(mesh, "censored depth", 0);
  assert(!vNumbering.isNull());
  const int oldNumVerticesNormal = vNumbering->size();

  int counterBegin = newNumCellsNormal + oldNumVerticesNormal;
  point_type curNewVertex = counterBegin;
  for(mesh_type::label_sequence::iterator c_iter = cellsNormal->begin(); c_iter != cellsNormalEnd; ++c_iter) {
    cV.clear();
    sieve->cone(*c_iter, cV);
    refiner.splitCell(*c_iter, cV.getPoints(), cV.getSize(), &curNewVertex);
  } // for
  for(mesh_type::label_sequence::iterator c_iter = cellsCensoredBegin; c_iter != cellsEnd; ++c_iter) {
    cV.clear();
    sieve->cone(*c_iter, cV);
    refiner.splitCell(*c_iter, cV.getPoints(), cV.getSize(), &curNewVertex);
  } // for
  const int newNumVerticesNormal = curNewVertex - counterBegin + oldNumVerticesNormal;
	
  // Count number of remaining vertices (other vertices).
  const int oldNumVerticesCensored = vertices->size() - oldNumVerticesNormal;
  counterBegin = newNumCellsNormal + newNumVerticesNormal + oldNumVerticesCensored;
  curNewVertex = counterBegin;
  for(mesh_type::label_sequence::iterator c_iter = cellsCensoredBegin; c_iter != cellsEnd; ++c_iter) {
    cV.clear();
    sieve->cone(*c_iter, cV);
    refiner.splitCellUncensored(*c_iter, cV.getPoints(), cV.getSize(), &curNewVertex);
  } // for
  const int newNumVerticesCensored = curNewVertex - counterBegin + oldNumVerticesCensored;

  _orderNewMesh->cellsNormal(0, newNumCellsNormal);
  _orderNewMesh->verticesNormal(newNumCellsNormal, newNumCellsNormal+newNumVerticesNormal);
  _orderNewMesh->verticesCensored(newNumCellsNormal+newNumVerticesNormal, newNumCellsNormal+newNumVerticesNormal+newNumVerticesCensored);
  _orderNewMesh->cellsCensored(newNumCellsNormal+newNumVerticesNormal+newNumVerticesCensored,
			       newNumCellsNormal+newNumVerticesNormal+newNumVerticesCensored+newNumCellsCensored);
  
  // Allocate chart for new sieve.
  const Obj<mesh_type::sieve_type>& newSieve = newMesh->getSieve();
  assert(!newSieve.isNull());
  newSieve->setChart(mesh_type::sieve_type::chart_type(0, _orderNewMesh->cellsCensored().max()));

  // Create new sieve with correct sizes for refined cells
  point_type curNewCell = _orderNewMesh->cellsNormal().min();
  interval_type::const_iterator oldCellsEnd = _orderOldMesh->cellsNormal().end();
  for (interval_type::const_iterator c_iter=_orderOldMesh->cellsNormal().begin(); c_iter != oldCellsEnd; ++c_iter) {
    // Set new cone and support sizes
    cV.clear();
    sieve->cone(*c_iter, cV);
    const point_type* cone = cV.getPoints();
    const int coneSize = cV.getSize();

    const point_type* newCells;
    int numNewCells = 0;
    refiner.getNewCells(&newCells, &numNewCells, *c_iter, cone, coneSize, *_orderOldMesh, *_orderNewMesh);

    for(int iCell=0; iCell < numNewCells; ++iCell, ++curNewCell) {
      newSieve->setConeSize(curNewCell, coneSize);
      for(int iVertex=0; iVertex < coneSize; ++iVertex) {
	newSieve->addSupportSize(newCells[iCell*coneSize+iVertex], 1);
      } // for
    } // for
  } // for
  // Reset current new cell value and loop over censored cells.
  curNewCell = _orderNewMesh->cellsCensored().min();
  oldCellsEnd = _orderOldMesh->cellsCensored().end();
  for (interval_type::const_iterator c_iter=_orderOldMesh->cellsCensored().begin(); c_iter != oldCellsEnd; ++c_iter) {
    // Set new cone and support sizes
    cV.clear();
    sieve->cone(*c_iter, cV);
    const point_type* cone = cV.getPoints();
    const int coneSize = cV.getSize();

    const point_type* newCells;
    int numNewCells = 0;
    refiner.getNewCells(&newCells, &numNewCells, *c_iter, cone, coneSize, *_orderOldMesh, *_orderNewMesh);

    for(int iCell=0; iCell < numNewCells; ++iCell, ++curNewCell) {
      newSieve->setConeSize(curNewCell, coneSize);
      for(int iVertex=0; iVertex < coneSize; ++iVertex) {
	newSieve->addSupportSize(newCells[iCell*coneSize+iVertex], 1);
      } // for
    } // for
  } // for
  newSieve->allocate();

  // Create refined cells in new sieve.
  curNewCell = _orderNewMesh->cellsNormal().min();
  oldCellsEnd = _orderOldMesh->cellsNormal().end();
  for (interval_type::const_iterator c_iter=_orderOldMesh->cellsNormal().begin(); c_iter != oldCellsEnd; ++c_iter) {
    cV.clear();
    sieve->cone(*c_iter, cV);
    const point_type *cone = cV.getPoints();
    const int coneSize = cV.getSize();
    
    const point_type* newCells;
    int numNewCells = 0;
    refiner.getNewCells(&newCells, &numNewCells, *c_iter, cone, coneSize, *_orderOldMesh, *_orderNewMesh);
    
    for(int iCell=0; iCell < numNewCells; ++iCell, ++curNewCell) {
      newSieve->setCone(&newCells[iCell*coneSize], curNewCell);
    } // for
  } // for
  curNewCell = _orderNewMesh->cellsCensored().min();
  oldCellsEnd = _orderOldMesh->cellsCensored().end();
  for (interval_type::const_iterator c_iter=_orderOldMesh->cellsCensored().begin(); c_iter != oldCellsEnd; ++c_iter) {
    cV.clear();
    sieve->cone(*c_iter, cV);
    const point_type *cone = cV.getPoints();
    const int coneSize = cV.getSize();
    
    const point_type* newCells;
    int numNewCells = 0;
    refiner.getNewCells(&newCells, &numNewCells, *c_iter, cone, coneSize, *_orderOldMesh, *_orderNewMesh);
    
    for(int iCell=0; iCell < numNewCells; ++iCell, ++curNewCell) {
      newSieve->setCone(&newCells[iCell*coneSize], curNewCell);
    } // for
  } // for
  newSieve->symmetrize();

  // Set coordinates in refined mesh.
  const Obj<mesh_type::real_section_type>& coordinates = mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const Obj<mesh_type::real_section_type>& newCoordinates = newMesh->getRealSection("coordinates");
  assert(!newCoordinates.isNull());

  const mesh_type::label_sequence::const_iterator verticesEnd = vertices->end();
  assert(vertices->size() > 0);
  const int spaceDim = coordinates->getFiberDimension(*vertices->begin());
  assert(spaceDim > 0);
  newCoordinates->setChart(mesh_type::sieve_type::chart_type(_orderNewMesh->verticesNormal().min(), _orderNewMesh->verticesCensored().max()));

  const interval_type::const_iterator newVerticesEnd = _orderNewMesh->verticesCensored().end();
  for (interval_type::const_iterator v_iter=_orderNewMesh->verticesNormal().begin(); v_iter != newVerticesEnd; ++v_iter) {
    newCoordinates->setFiberDimension(*v_iter, spaceDim);
  } // for
  newCoordinates->allocatePoint();
      
  interval_type::const_iterator oldVerticesEnd = _orderOldMesh->verticesNormal().end();
  for (interval_type::const_iterator vOld_iter=_orderOldMesh->verticesNormal().begin(), vNew_iter=_orderNewMesh->verticesNormal().begin(); vOld_iter != oldVerticesEnd; ++vOld_iter, ++vNew_iter) {
    //std::cout << "Copy coordinates from old vertex " << *vOld_iter << " to new vertex " << *vNew_iter << std::endl;
    newCoordinates->updatePoint(*vNew_iter, coordinates->restrictPoint(*vOld_iter));
  } // for
  oldVerticesEnd = _orderOldMesh->verticesCensored().end();
  for (interval_type::const_iterator vOld_iter=_orderOldMesh->verticesCensored().begin(), vNew_iter=_orderNewMesh->verticesCensored().begin(); vOld_iter != oldVerticesEnd; ++vOld_iter, ++vNew_iter) {
    //std::cout << "Copy coordinates from old vertex " << *vOld_iter << " to new vertex " << *vNew_iter << std::endl;
    newCoordinates->updatePoint(*vNew_iter, coordinates->restrictPoint(*vOld_iter));
  } // for

  refiner.setCoordsNewVertices(newCoordinates, coordinates);

  _stratify(newMesh);
  _calcNewOverlap(newMesh, mesh, refiner);
  _createIntSections(newMesh, mesh, refiner);
  _createLabels(newMesh, mesh, refiner);

  // Create sensored depth
  const ALE::Obj<Mesh::label_type>& censoredLabel = newMesh->createLabel("censored depth");
  assert(!censoredLabel.isNull());

  mesh_type::DepthVisitor depthVisitor(*newSieve, _orderNewMesh->verticesCensored().min(), *censoredLabel);

  newSieve->roots(depthVisitor);
  while(depthVisitor.isModified()) {
    // FIX: Avoid the copy here somehow by fixing the traversal
    std::vector<mesh_type::point_type> modifiedPoints(depthVisitor.getModifiedPoints().begin(), depthVisitor.getModifiedPoints().end());
    
    depthVisitor.clear();
    newSieve->support(modifiedPoints, depthVisitor);
  } // while
} // _refineCensored

// ----------------------------------------------------------------------
// Stratify mesh.
template<typename cellrefiner_type>
void
ALE::MeshRefiner<cellrefiner_type>::_stratify(const Obj<mesh_type>& mesh)
{ // _stratify
  typedef Interval<point_type> interval_type;

  assert(_orderNewMesh);
  assert(_orderOldMesh);
	 
  // Fast stratification
  const Obj<mesh_type::label_type>& height = mesh->createLabel("height");
  const Obj<mesh_type::label_type>& depth  = mesh->createLabel("depth");

  // Set height/depth of cells
  interval_type::const_iterator cellsEnd = _orderNewMesh->cellsNormal().end();
  for (interval_type::const_iterator c_iter = _orderNewMesh->cellsNormal().begin(); c_iter != cellsEnd; ++c_iter) {
    height->setCone(0, *c_iter);
    depth->setCone(1, *c_iter);
  } // for
  cellsEnd = _orderNewMesh->cellsCensored().end();
  for (interval_type::const_iterator c_iter = _orderNewMesh->cellsCensored().begin(); c_iter != cellsEnd; ++c_iter) {
    height->setCone(0, *c_iter);
    depth->setCone(1, *c_iter);
  } // for

  // Set height/depth of vertices
  interval_type::const_iterator verticesEnd = _orderNewMesh->verticesNormal().end();
  for (interval_type::const_iterator v_iter = _orderNewMesh->verticesNormal().begin(); v_iter != verticesEnd; ++v_iter) {
    height->setCone(1, *v_iter);
    depth->setCone(0, *v_iter);
  } // for
  verticesEnd = _orderNewMesh->verticesCensored().end();
  for (interval_type::const_iterator v_iter = _orderNewMesh->verticesCensored().begin(); v_iter != verticesEnd; ++v_iter) {
    height->setCone(1, *v_iter);
    depth->setCone(0, *v_iter);
  } // for
  
  mesh->setHeight(1);
  mesh->setDepth(1);
} // _stratify

// ----------------------------------------------------------------------
// Create integer sections in new mesh.
template<typename cellrefiner_type>
void
ALE::MeshRefiner<cellrefiner_type>::_createIntSections(const Obj<mesh_type>& newMesh,
						  const Obj<mesh_type>& mesh,
						  cellrefiner_type& refiner)
{ // _createIntSections
  assert(!newMesh.isNull());
  assert(!mesh.isNull());

  const ALE::Obj<std::set<std::string> >& sectionNames =
    mesh->getIntSections();  
  const std::set<std::string>::const_iterator namesBegin = 
    sectionNames->begin();
  const std::set<std::string>::const_iterator namesEnd = 
    sectionNames->end();
  for (std::set<std::string>::const_iterator name=namesBegin;
      name != namesEnd;
      ++name) {
    const ALE::Obj<mesh_type::int_section_type>& group = mesh->getIntSection(*name);
    const ALE::Obj<mesh_type::int_section_type>& newGroup = newMesh->getIntSection(*name);
    const mesh_type::int_section_type::chart_type& chart = group->getChart();

    // :WARNING: Only implemented for integer sections containing vertices.

    const point_type oldVerticesStart = _orderOldMesh->verticesNormal().min();
    const point_type oldVerticesCensoredStart = _orderOldMesh->verticesCensored().min();
    const point_type oldVerticesEnd = _orderOldMesh->verticesCensored().max();
    
    const point_type newVerticesStart = _orderNewMesh->verticesNormal().min();
    const point_type newVerticesCensoredStart = _orderNewMesh->verticesCensored().min();
    const point_type newVerticesEnd = _orderNewMesh->verticesCensored().max();
    
    newGroup->setChart(mesh_type::int_section_type::chart_type(newVerticesStart, newVerticesEnd));
    const mesh_type::int_section_type::chart_type& newChart = newGroup->getChart();

    //std::cout << "VERTICES start: " << newVerticesStart << ", end: " << newVerticesEnd << std::endl;

    const point_type chartMax = chart.max();
    for (point_type pOld=chart.min(); pOld < chartMax; ++pOld) {
      if (group->getFiberDimension(pOld)) {
	if (_orderOldMesh->verticesNormal().hasPoint(pOld)) {
	  const point_type pNew = newVerticesStart + pOld - oldVerticesStart;
	  newGroup->setFiberDimension(pNew, 1);
	} else if (_orderOldMesh->verticesCensored().hasPoint(pOld)) {
	  const point_type pNew = newVerticesCensoredStart + pOld - oldVerticesCensoredStart;
	  newGroup->setFiberDimension(pNew, 1);
	} else {
	  throw ALE::Exception("Creating integer sections during refinement containing entities other than vertices not implemented.");
	} // if/else
      } // if
    } // for
    refiner.groupAddNewVertices(newGroup, group);

    newGroup->allocatePoint();
    for (point_type pOld=chart.min(); pOld < chartMax; ++pOld) {
      if (group->getFiberDimension(pOld)) {
	if (_orderOldMesh->verticesNormal().hasPoint(pOld)) {
	  const point_type pNew = newVerticesStart + pOld - oldVerticesStart;
	  newGroup->updatePoint(pNew, group->restrictPoint(pOld));
	} else if (_orderOldMesh->verticesCensored().hasPoint(pOld)) {
	  const point_type pNew = newVerticesCensoredStart + pOld - oldVerticesCensoredStart;
	  newGroup->updatePoint(pNew, group->restrictPoint(pOld));
	} else {
	  throw ALE::Exception("Creating integer sections during refinement containing entities other than vertices not implemented.");
	} // if/else
      } // if
    } // for
    refiner.groupSetNewVertices(newGroup, group);
  } // for
} // _createIntSections

// ----------------------------------------------------------------------
// Create labels in new mesh.
template<typename cellrefiner_type>
void
ALE::MeshRefiner<cellrefiner_type>::_createLabels(const Obj<mesh_type>& newMesh,
					     const Obj<mesh_type>& mesh,
					     cellrefiner_type& refiner)
{ // _createLabels
  assert(!newMesh.isNull());
  assert(!mesh.isNull());
  assert(_orderOldMesh);
  assert(_orderNewMesh);

  typedef ALE::Interval<point_type> interval_type;

  const mesh_type::labels_type labels = mesh->getLabels();
  const mesh_type::labels_type::const_iterator labelsEnd = labels.end();
  for (mesh_type::labels_type::const_iterator l_iter=labels.begin(); l_iter != labelsEnd; ++l_iter) {
    // Handle censored depth separately.
    if ("censored depth" == l_iter->first || "depth" == l_iter->first || "height" == l_iter->first) {
      continue;
    } // if

    const Obj<mesh_type::label_type>& oldLabel = l_iter->second;
    assert(!oldLabel.isNull());
    const Obj<mesh_type::label_type>& newLabel = newMesh->createLabel(l_iter->first);
    assert(!newLabel.isNull());

    const int defaultValue = -999;

    // Update cells
    // Normal cells
    interval_type::const_iterator pointsEnd = _orderOldMesh->cellsNormal().end();
    for (interval_type::const_iterator p_iter = _orderOldMesh->cellsNormal().begin(); p_iter != pointsEnd; ++p_iter) {
      const int value = mesh->getValue(oldLabel, *p_iter, defaultValue);
      if (defaultValue == value)
	continue;

      const int numNewCellsPerCell = refiner.numNewCells(*p_iter);
      mesh_type::point_type pNew = _orderNewMesh->cellsNormal().min() + (*p_iter - _orderOldMesh->cellsNormal().min())*numNewCellsPerCell;
      for(int i=0; i < numNewCellsPerCell; ++i, ++pNew)
	newMesh->setValue(newLabel, pNew, value);
    } // for
    
    // Censored cells
    pointsEnd = _orderOldMesh->cellsCensored().end();
    for (interval_type::const_iterator p_iter = _orderOldMesh->cellsCensored().begin(); p_iter != pointsEnd; ++p_iter) {
      const int value = mesh->getValue(oldLabel, *p_iter, defaultValue);
      if (defaultValue == value)
	continue;
      
      const int numNewCellsPerCell = refiner.numNewCells(*p_iter);
      mesh_type::point_type pNew = _orderNewMesh->cellsCensored().min() + (*p_iter - _orderOldMesh->cellsCensored().min())*numNewCellsPerCell;
      for(int i=0; i < numNewCellsPerCell; ++i, ++pNew)
	newMesh->setValue(newLabel, pNew, value);
    } // for
    
    // Normal vertices
    pointsEnd = _orderOldMesh->verticesNormal().end();
    for (interval_type::const_iterator p_iter = _orderOldMesh->verticesNormal().begin(); p_iter != pointsEnd; ++p_iter) {
      const int value = mesh->getValue(oldLabel, *p_iter, defaultValue);
      if (defaultValue == value)
	continue;
      
      const mesh_type::point_type pNew = _orderNewMesh->verticesNormal().min() + (*p_iter - _orderOldMesh->verticesNormal().min());
      newMesh->setValue(newLabel, pNew, value);
    } // for

    // Censored vertices
    pointsEnd = _orderOldMesh->verticesCensored().end();
    for (interval_type::const_iterator p_iter = _orderOldMesh->verticesCensored().begin(); p_iter != pointsEnd; ++p_iter) {
      const int value = mesh->getValue(oldLabel, *p_iter, defaultValue);
      if (defaultValue == value)
	continue;
      
      const mesh_type::point_type pNew = _orderNewMesh->verticesCensored().min() + (*p_iter - _orderOldMesh->verticesCensored().min());
      newMesh->setValue(newLabel, pNew, value);
    } // for

    refiner.labelAddNewVertices(newMesh, mesh, l_iter->first.c_str());
  } // for
} // _createLabels

// ----------------------------------------------------------------------
// Calculate new overlap.
template<typename cellrefiner_type>
void
ALE::MeshRefiner<cellrefiner_type>::_calcNewOverlap(const Obj<mesh_type>& newMesh,
						    const Obj<mesh_type>& mesh,
						    cellrefiner_type& refiner)
{ // _calcNewOverlap
  assert(!newMesh.isNull());
  assert(!mesh.isNull());

  // Exchange new boundary vertices
  //   We can convert endpoints, and then just match to new vertex on this side
  //   1) Create the overlap of edges which are vertex pairs (do not need for interpolated meshes)
  //   2) Create a section of overlap edge --> new vertex (this will generalize to other split points in interpolated meshes)
  //   3) Copy across new overlap
  //   4) Fuse matches new vertex pairs and inserts them into the old overlap

  // Create the parallel overlap

  // Get offsets for points across processors and add points in overlap from original mesh to new overlap
  int* oldVerticesStartP = (mesh->commSize() > 0) ? new int[mesh->commSize()] : 0;
  int* newVerticesStartP = (newMesh->commSize() > 0) ? new int[newMesh->commSize()] : 0;
  int ierr = 0;
  
  const int oldVerticesStart = _orderOldMesh->verticesNormal().min();
  const int newVerticesStart = _orderNewMesh->verticesNormal().min();

  ierr = MPI_Allgather((void *) &oldVerticesStart, 1, MPI_INT, oldVerticesStartP, 1, MPI_INT, mesh->comm());CHKERRXX(ierr);
  ierr = MPI_Allgather((void *) &newVerticesStart, 1, MPI_INT, newVerticesStartP, 1, MPI_INT, newMesh->comm());CHKERRXX(ierr);
  Obj<mesh_type::send_overlap_type> newSendOverlap = newMesh->getSendOverlap();
  Obj<mesh_type::recv_overlap_type> newRecvOverlap = newMesh->getRecvOverlap();
  const Obj<mesh_type::send_overlap_type>& sendOverlap = mesh->getSendOverlap();
  const Obj<mesh_type::recv_overlap_type>& recvOverlap = mesh->getRecvOverlap();
  Obj<mesh_type::send_overlap_type::traits::capSequence> sendPoints  = sendOverlap->cap();
  const mesh_type::send_overlap_type::source_type localOffset = newVerticesStart - oldVerticesStart;
  
  for(mesh_type::send_overlap_type::traits::capSequence::iterator p_iter = sendPoints->begin(); p_iter != sendPoints->end(); ++p_iter) {
    const Obj<mesh_type::send_overlap_type::traits::supportSequence>& ranks      = sendOverlap->support(*p_iter);
    const mesh_type::send_overlap_type::source_type&                  localPoint = *p_iter;
    
    for(mesh_type::send_overlap_type::traits::supportSequence::iterator r_iter = ranks->begin(); r_iter != ranks->end(); ++r_iter) {
      const int                                   rank         = *r_iter;
      const mesh_type::send_overlap_type::source_type& remotePoint  = r_iter.color();
      const mesh_type::send_overlap_type::source_type  remoteOffset = newVerticesStartP[rank] - oldVerticesStartP[rank];
      
      newSendOverlap->addArrow(localPoint+localOffset, rank, remotePoint+remoteOffset);
    } // for
  } // for
  Obj<mesh_type::recv_overlap_type::traits::baseSequence> recvPoints = recvOverlap->base();
  
  for(mesh_type::recv_overlap_type::traits::baseSequence::iterator p_iter = recvPoints->begin(); p_iter != recvPoints->end(); ++p_iter) {
    const Obj<mesh_type::recv_overlap_type::traits::coneSequence>& ranks      = recvOverlap->cone(*p_iter);
    const mesh_type::recv_overlap_type::target_type&               localPoint = *p_iter;
    
    for(mesh_type::recv_overlap_type::traits::coneSequence::iterator r_iter = ranks->begin(); r_iter != ranks->end(); ++r_iter) {
      const int                                        rank         = *r_iter;
      const mesh_type::recv_overlap_type::target_type& remotePoint  = r_iter.color();
      const mesh_type::recv_overlap_type::target_type  remoteOffset = newVerticesStartP[rank] - oldVerticesStartP[rank];
      
      newRecvOverlap->addArrow(rank, localPoint+localOffset, remotePoint+remoteOffset);
    } // for
  } // for
  newMesh->setCalculatedOverlap(true);
  delete [] oldVerticesStartP; oldVerticesStartP = PETSC_NULL;
  delete [] newVerticesStartP; newVerticesStartP = PETSC_NULL;

  refiner.overlapAddNewVertices(newMesh, *_orderNewMesh, mesh, *_orderOldMesh);
} // _calcNewOverlap


// End of file 
