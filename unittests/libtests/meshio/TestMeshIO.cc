// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestMeshIO.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/meshio/MeshIO.hh" // USES MeshIO

#include "pylith/utils/array.hh" // USES int_array

#include "data/MeshData.hh"

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::logic_error

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Get simple mesh for testing I/O.
pylith::topology::Mesh*
pylith::meshio::TestMeshIO::_createMesh(const MeshData& data)
{ // _createMesh
  PYLITH_METHOD_BEGIN;

  // buildTopology() requires zero based index
  CPPUNIT_ASSERT(true == data.useIndexZero);

  CPPUNIT_ASSERT(data.vertices);
  CPPUNIT_ASSERT(data.cells);
  CPPUNIT_ASSERT(data.materialIds);
  if (data.numGroups > 0) {
    CPPUNIT_ASSERT(data.groups);
    CPPUNIT_ASSERT(data.groupSizes);
    CPPUNIT_ASSERT(data.groupNames);
    CPPUNIT_ASSERT(data.groupTypes);
  } // if

  topology::Mesh* mesh = new topology::Mesh(data.cellDim);CPPUNIT_ASSERT(mesh);
  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();CPPUNIT_ASSERT(!sieveMesh.isNull());
  ALE::Obj<SieveMesh::sieve_type> sieve = new SieveMesh::sieve_type(mesh->comm());CPPUNIT_ASSERT(!sieve.isNull());

  // Cells and vertices
  const bool interpolate = false;
  ALE::Obj<SieveFlexMesh::sieve_type> s = new SieveFlexMesh::sieve_type(sieve->comm(), sieve->debug());CPPUNIT_ASSERT(s);
  ALE::SieveBuilder<SieveFlexMesh>::buildTopology(s, data.cellDim, data.numCells, const_cast<int*>(data.cells), data.numVertices, interpolate, data.numCorners);
  std::map<SieveFlexMesh::point_type,SieveFlexMesh::point_type> renumbering;
  ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
  sieveMesh->setSieve(sieve);
  sieveMesh->stratify();
  ALE::SieveBuilder<SieveMesh>::buildCoordinates(sieveMesh, data.spaceDim, data.vertices);

  PetscDM dmMesh = NULL;
  PetscBool interpolateMesh = interpolate ? PETSC_TRUE : PETSC_FALSE;
  PetscErrorCode err;

  err = DMPlexCreateFromCellList(mesh->comm(), data.cellDim, data.numCells, data.numVertices, data.numCorners, interpolateMesh, data.cells, data.spaceDim, data.vertices, &dmMesh);CHECK_PETSC_ERROR(err);
  mesh->setDMMesh(dmMesh);

  // Material ids
  const ALE::Obj<SieveMesh::label_sequence>& cells = sieveMesh->heightStratum(0);
  CPPUNIT_ASSERT(!cells.isNull());
  const ALE::Obj<SieveMesh::label_type>& labelMaterials = sieveMesh->createLabel("material-id");
  CPPUNIT_ASSERT(!labelMaterials.isNull());
  int i = 0;
  for(SieveMesh::label_sequence::iterator e_iter=cells->begin(); e_iter != cells->end(); ++e_iter)
    sieveMesh->setValue(labelMaterials, *e_iter, data.materialIds[i++]);

  PetscInt cStart, cEnd;
  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    err = DMPlexSetLabelValue(dmMesh, "material-id", c, data.materialIds[c-cStart]);CHECK_PETSC_ERROR(err);
  } // for

  // Groups
  for (int iGroup=0, index=0; iGroup < data.numGroups; ++iGroup) {
    const ALE::Obj<SieveMesh::int_section_type>& groupField = sieveMesh->getIntSection(data.groupNames[iGroup]);
    CPPUNIT_ASSERT(!groupField.isNull());
    groupField->setChart(sieveMesh->getSieve()->getChart());

    err = DMPlexCreateLabel(dmMesh, data.groupNames[iGroup]);CHECK_PETSC_ERROR(err);

    MeshIO::GroupPtType type;
    const int numPoints = data.groupSizes[iGroup];
    if (0 == strcasecmp("cell", data.groupTypes[iGroup])) {
      type = MeshIO::CELL;
      for(int i=0; i < numPoints; ++i, ++index) {
        groupField->setFiberDimension(data.groups[index], 1);
        err = DMPlexSetLabelValue(dmMesh, data.groupNames[iGroup], data.groups[index], 1);CHECK_PETSC_ERROR(err);
      } // for
    } else if (0 == strcasecmp("vertex", data.groupTypes[iGroup])) {
      type = MeshIO::VERTEX;
      const int numCells = sieveMesh->heightStratum(0)->size();
      for(int i=0; i < numPoints; ++i, ++index) {
        groupField->setFiberDimension(data.groups[index]+numCells, 1);
        err = DMPlexSetLabelValue(dmMesh, data.groupNames[iGroup], data.groups[index]+numCells, 1);CHECK_PETSC_ERROR(err);
      } // for
    } else
      throw std::logic_error("Could not parse group type.");
    sieveMesh->allocate(groupField);
  } // for
 
    // Set coordinate system
  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(data.spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);

  PYLITH_METHOD_RETURN(mesh);
} // _createMesh

// ----------------------------------------------------------------------
// Check values in mesh against data.
void
pylith::meshio::TestMeshIO::_checkVals(const topology::Mesh& mesh,
				       const MeshData& data)
{ // _checkVals
  PYLITH_METHOD_BEGIN;

  typedef SieveMesh::int_section_type::chart_type chart_type;

  // Check mesh dimension
  CPPUNIT_ASSERT_EQUAL(data.cellDim, mesh.dimension());

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();CPPUNIT_ASSERT(!sieveMesh.isNull());

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);

  // Check vertices
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  const ALE::Obj<RealSection>& coordsField =
    sieveMesh->getRealSection("coordinates");
  const int numVertices = vertices->size();
  CPPUNIT_ASSERT(!vertices.isNull());
  CPPUNIT_ASSERT(!coordsField.isNull());
  CPPUNIT_ASSERT_EQUAL(data.numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(data.spaceDim, coordsField->getFiberDimension(*vertices->begin()));
  int i = 0;
  const int spaceDim = data.spaceDim;
  for(SieveMesh::label_sequence::iterator v_iter = vertices->begin(); v_iter != vertices->end(); ++v_iter) {
    const PylithScalar* vertexCoords = coordsField->restrictPoint(*v_iter);CPPUNIT_ASSERT(vertexCoords);
    const PylithScalar tolerance = 1.0e-06;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      if (data.vertices[i] < 1.0) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(data.vertices[i++], vertexCoords[iDim],
				     tolerance);
      } else {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vertexCoords[iDim]/data.vertices[i++],
				     tolerance);
      } // if/else
  } // for
  PetscSection coordSection = NULL;
  PetscVec coordVec = NULL;
  PetscScalar *coords = NULL;
  PetscInt vStart, vEnd;
  PetscErrorCode err;

  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);CPPUNIT_ASSERT(coordSection);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);CPPUNIT_ASSERT(coordVec);
  err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off;

    err = PetscSectionGetDof(coordSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(coordSection, v, &off);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dof);
    const PylithScalar tolerance = 1.0e-06;
    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PetscInt eoff = (v-vStart)*spaceDim;
      if (data.vertices[eoff+d] < 1.0) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(data.vertices[eoff+d], coords[off+d], tolerance);
      } else {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, coords[off+d]/data.vertices[eoff+d], tolerance);
      } // if/else
    } // for
  } // for
  err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);

  // check cells
  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  const ALE::Obj<SieveMesh::label_sequence>& cells = sieveMesh->heightStratum(0);

  const int numCells = cells->size();
  CPPUNIT_ASSERT_EQUAL(data.numCells, numCells);
  const int numCorners = sieveMesh->getNumCellCorners();
  CPPUNIT_ASSERT_EQUAL(data.numCorners, numCorners);

  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> pV(sieve->getMaxConeSize());
  const int offset = numCells;
  i = 0;
  for(SieveMesh::label_sequence::iterator e_iter = cells->begin();
      e_iter != cells->end();
      ++e_iter) {
    sieve->cone(*e_iter, pV);
    const SieveMesh::point_type *cone = pV.getPoints();
    for(int p = 0; p < pV.getSize(); ++p, ++i) {
      CPPUNIT_ASSERT_EQUAL(data.cells[i], cone[p]-offset);
    } // for
    pV.clear();
  } // for
  PetscInt cStart, cEnd;

  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(data.numCells, cEnd-cStart);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    const PetscInt *cone;
    PetscInt        coneSize;

    err = DMPlexGetConeSize(dmMesh, c, &coneSize);CHECK_PETSC_ERROR(err);
    err = DMPlexGetCone(dmMesh, c, &cone);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(data.numCorners, coneSize);
    for(PetscInt p = 0; p < coneSize; ++p) {
      CPPUNIT_ASSERT_EQUAL(data.cells[(c-cStart)*coneSize+p], cone[p]-offset);
    } // for
  } // for

  // check materials
  const ALE::Obj<SieveMesh::label_type>& labelMaterials = sieveMesh->getLabel("material-id");
  const int idDefault = -999;
  const int size = numCells;
  int_array materialIds(size);
  i = 0;
  for(SieveMesh::label_sequence::iterator e_iter = cells->begin();
      e_iter != cells->end();
      ++e_iter)
    materialIds[i++] = sieveMesh->getValue(labelMaterials, *e_iter, idDefault);
  
  for (int iCell=0; iCell < numCells; ++iCell)
    CPPUNIT_ASSERT_EQUAL(data.materialIds[iCell], materialIds[iCell]);

  for(PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt matId;

    err = DMPlexGetLabelValue(dmMesh, "material-id", c, &matId);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(data.materialIds[c-cStart], matId);
  } // for

  // Check groups
#if 0 // SIEVE STUFF
  const ALE::Obj<std::set<std::string> >& groupNames = 
    sieveMesh->getIntSections();
  if (data.numGroups > 0) {
    CPPUNIT_ASSERT(!groupNames.isNull());
    CPPUNIT_ASSERT_EQUAL(data.numGroups, int(groupNames->size()));
  } // if
  int iGroup = 0;
  int index = 0;
  for (std::set<std::string>::const_iterator name=groupNames->begin();
       name != groupNames->end();
       ++name, ++iGroup) {
    const ALE::Obj<SieveMesh::int_section_type>& groupField = 
      sieveMesh->getIntSection(*name);
    CPPUNIT_ASSERT(!groupField.isNull());
    const chart_type& chart = groupField->getChart();
    SieveMesh::point_type firstPoint;
    for(chart_type::const_iterator c_iter = chart.begin();
	c_iter != chart.end();
	++c_iter) {
      if (groupField->getFiberDimension(*c_iter)) {
	firstPoint = *c_iter;
	break;
      } // if
    } // for
    std::string groupType = 
      (sieveMesh->height(firstPoint) == 0) ? "cell" : "vertex";
    const int numPoints = groupField->size();
    int_array points(numPoints);
    int i = 0;
    const int offset = ("vertex" == groupType) ? numCells : 0;
    for(chart_type::const_iterator c_iter = chart.begin();
	c_iter != chart.end();
	++c_iter)
      if (groupField->getFiberDimension(*c_iter))
	points[i++] = *c_iter - offset;
    
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupNames[iGroup]), *name);
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupTypes[iGroup]), groupType);
    CPPUNIT_ASSERT_EQUAL(data.groupSizes[iGroup], numPoints);
    for (int i=0; i < numPoints; ++i)
      CPPUNIT_ASSERT_EQUAL(data.groups[index++], points[i]);
  } // for
#endif

  PetscInt numLabels, pStart, pEnd;

  err = DMPlexGetChart(dmMesh, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetNumLabels(dmMesh, &numLabels);CHECK_PETSC_ERROR(err);
  numLabels -= 2; /* Remove depth and material labels */
  CPPUNIT_ASSERT_EQUAL(data.numGroups, numLabels);
  PetscInt iGroup = 0;
  PetscInt index  = 0;
  for(PetscInt l = numLabels-1; l >= 0; --l, ++iGroup) {
    const char *name;
    PetscInt    firstPoint;

    err = DMPlexGetLabelName(dmMesh, l, &name);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupNames[iGroup]), std::string(name));
    for(PetscInt p = pStart; p < pEnd; ++p) {
      PetscInt val;

      err = DMPlexGetLabelValue(dmMesh, name, p, &val);CHECK_PETSC_ERROR(err);
      if (val >= 0) {
        firstPoint = p;
        break;
      } // if
    } // for
    std::string groupType = (firstPoint >= cStart && firstPoint < cEnd) ? "cell" : "vertex";
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupTypes[iGroup]), groupType);
    PetscInt numPoints;
    err = DMPlexGetStratumSize(dmMesh, name, 1, &numPoints);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(data.groupSizes[iGroup], numPoints);
    PetscIS pointIS = NULL;
    const PetscInt *points = NULL;
    const PetscInt offset = ("vertex" == groupType) ? numCells : 0;
    err = DMPlexGetStratumIS(dmMesh, name, 1, &pointIS);CHECK_PETSC_ERROR(err);
    err = ISGetIndices(pointIS, &points);CHECK_PETSC_ERROR(err);
    for(PetscInt p = 0; p < numPoints; ++p) {
      CPPUNIT_ASSERT_EQUAL(data.groups[index++], points[p]-offset);
    } // for
    err = ISRestoreIndices(pointIS, &points);CHECK_PETSC_ERROR(err);
    err = ISDestroy(&pointIS);CHECK_PETSC_ERROR(err);
  } // for

  PYLITH_METHOD_END;
} // _checkVals

// ----------------------------------------------------------------------
// Test debug()
void
pylith::meshio::TestMeshIO::_testDebug(MeshIO& iohandler)
{ // _testDebug
  PYLITH_METHOD_BEGIN;

  bool debug = false;
  iohandler.debug(debug);
  CPPUNIT_ASSERT_EQUAL(debug, iohandler.debug());
  
  debug = true;
  iohandler.debug(debug);
  CPPUNIT_ASSERT_EQUAL(debug, iohandler.debug());  

  PYLITH_METHOD_END;
} // _testDebug

// ----------------------------------------------------------------------
// Test interpolate()
void
pylith::meshio::TestMeshIO::_testInterpolate(MeshIO& iohandler)
{ // _testInterpolate
  PYLITH_METHOD_BEGIN;

  bool interpolate = false;
  iohandler.interpolate(interpolate);
  CPPUNIT_ASSERT_EQUAL(interpolate, iohandler.interpolate());
  
  interpolate = true;
  iohandler.interpolate(interpolate);
  CPPUNIT_ASSERT_EQUAL(interpolate, iohandler.interpolate());  

  PYLITH_METHOD_END;
} // _testInterpolate


// End of file 
