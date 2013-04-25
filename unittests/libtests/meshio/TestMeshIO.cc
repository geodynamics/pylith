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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestMeshIO.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/meshio/MeshIO.hh" // USES MeshIO

#include "pylith/utils/array.hh" // USES int_array

#include "data/MeshData.hh"

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::logic_error

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestMeshIO::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _mesh = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestMeshIO::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _mesh; _mesh = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Get simple mesh for testing I/O.
void
pylith::meshio::TestMeshIO::_createMesh(const MeshData& data)
{ // _createMesh
  PYLITH_METHOD_BEGIN;

  // buildTopology() requires zero based index
  CPPUNIT_ASSERT_EQUAL(true, data.useIndexZero);

  CPPUNIT_ASSERT(data.vertices);
  CPPUNIT_ASSERT(data.cells);
  CPPUNIT_ASSERT(data.materialIds);
  if (data.numGroups > 0) {
    CPPUNIT_ASSERT(data.groups);
    CPPUNIT_ASSERT(data.groupSizes);
    CPPUNIT_ASSERT(data.groupNames);
    CPPUNIT_ASSERT(data.groupTypes);
  } // if

  delete _mesh; _mesh = new topology::Mesh(data.cellDim);CPPUNIT_ASSERT(_mesh);

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();CPPUNIT_ASSERT(!sieveMesh.isNull());
  ALE::Obj<SieveMesh::sieve_type> sieve = new SieveMesh::sieve_type(_mesh->comm());CPPUNIT_ASSERT(!sieve.isNull());

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

  err = DMPlexCreateFromCellList(_mesh->comm(), data.cellDim, data.numCells, data.numVertices, data.numCorners, interpolateMesh, data.cells, data.spaceDim, data.vertices, &dmMesh);PYLITH_CHECK_ERROR(err);
  _mesh->setDMMesh(dmMesh);

  // Material ids
  const ALE::Obj<SieveMesh::label_sequence>& cells = sieveMesh->heightStratum(0);
  CPPUNIT_ASSERT(!cells.isNull());
  const ALE::Obj<SieveMesh::label_type>& labelMaterials = sieveMesh->createLabel("material-id");
  CPPUNIT_ASSERT(!labelMaterials.isNull());
  int i = 0;
  for(SieveMesh::label_sequence::iterator e_iter=cells->begin(); e_iter != cells->end(); ++e_iter)
    sieveMesh->setValue(labelMaterials, *e_iter, data.materialIds[i++]);

  PetscInt cStart, cEnd;
  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    err = DMPlexSetLabelValue(dmMesh, "material-id", c, data.materialIds[c-cStart]);PYLITH_CHECK_ERROR(err);
  } // for

  // Groups
  for (int iGroup=0, index=0; iGroup < data.numGroups; ++iGroup) {
    const ALE::Obj<SieveMesh::int_section_type>& groupField = sieveMesh->getIntSection(data.groupNames[iGroup]);
    CPPUNIT_ASSERT(!groupField.isNull());
    groupField->setChart(sieveMesh->getSieve()->getChart());

    err = DMPlexCreateLabel(dmMesh, data.groupNames[iGroup]);PYLITH_CHECK_ERROR(err);

    MeshIO::GroupPtType type;
    const int numPoints = data.groupSizes[iGroup];
    if (0 == strcasecmp("cell", data.groupTypes[iGroup])) {
      type = MeshIO::CELL;
      for(int i=0; i < numPoints; ++i, ++index) {
        groupField->setFiberDimension(data.groups[index], 1);
        err = DMPlexSetLabelValue(dmMesh, data.groupNames[iGroup], data.groups[index], 1);PYLITH_CHECK_ERROR(err);
      } // for
    } else if (0 == strcasecmp("vertex", data.groupTypes[iGroup])) {
      type = MeshIO::VERTEX;
      const int numCells = sieveMesh->heightStratum(0)->size();
      for(int i=0; i < numPoints; ++i, ++index) {
        groupField->setFiberDimension(data.groups[index]+numCells, 1);
        err = DMPlexSetLabelValue(dmMesh, data.groupNames[iGroup], data.groups[index]+numCells, 1);PYLITH_CHECK_ERROR(err);
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
  _mesh->coordsys(&cs);
  _mesh->nondimensionalize(normalizer);

  PYLITH_METHOD_END;
} // _createMesh

// ----------------------------------------------------------------------
// Check values in mesh against data.
void
pylith::meshio::TestMeshIO::_checkVals(const MeshData& data)
{ // _checkVals
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);

  // Check mesh dimension
  CPPUNIT_ASSERT_EQUAL(data.cellDim, _mesh->dimension());
  const int spaceDim = data.spaceDim;

  PetscDM dmMesh = _mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);

  // Check vertices
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  CPPUNIT_ASSERT_EQUAL(data.numVertices, verticesStratum.size());

  topology::CoordsVisitor coordsVisitor(dmMesh);
  const PetscScalar* coordsArray = coordsVisitor.localArray();
  const PylithScalar tolerance = 1.0e-06;
  for (PetscInt v = vStart, index = 0; v < vEnd; ++v) {
    const PetscInt off = coordsVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, coordsVisitor.sectionDof(v));

    for (int iDim=0; iDim < spaceDim; ++iDim, ++index) {
      if (data.vertices[index] < 1.0) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(data.vertices[index], coordsArray[off+iDim], tolerance);
      } else {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, coordsArray[off+iDim]/data.vertices[index], tolerance);
      } // if/else
    } // for
  } // for

  // Check cells
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  const PetscInt numCells = cellsStratum.size();

  CPPUNIT_ASSERT_EQUAL(data.numCells, numCells);
  const int offset = numCells;
  const PetscInt* cone = NULL;
  PetscInt coneSize = 0;
  PetscErrorCode err = 0;
  for(PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    err = DMPlexGetConeSize(dmMesh, c, &coneSize);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetCone(dmMesh, c, &cone);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(data.numCorners, coneSize);
    for(PetscInt p = 0; p < coneSize; ++p, ++index) {
      CPPUNIT_ASSERT_EQUAL(data.cells[index], cone[p]-offset);
    } // for
  } // for

  // check materials
  PetscInt matId = 0;
  for(PetscInt c = cStart; c < cEnd; ++c) {
    err = DMPlexGetLabelValue(dmMesh, "material-id", c, &matId);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(data.materialIds[c-cStart], matId);
  } // for

  // Check groups
  PetscInt numGroups, pStart, pEnd;
  err = DMPlexGetChart(dmMesh, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetNumLabels(dmMesh, &numGroups);PYLITH_CHECK_ERROR(err);
  numGroups -= 2; // Remove depth and material labels.
  CPPUNIT_ASSERT_EQUAL(data.numGroups, numGroups);
  PetscInt index  = 0;
  for(PetscInt iGroup = 0, iLabel = numGroups-1; iGroup < numGroups; ++iGroup, --iLabel) {
    const char *name = NULL;
    PetscInt firstPoint = 0;

    err = DMPlexGetLabelName(dmMesh, iLabel, &name);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupNames[iGroup]), std::string(name));
    for(PetscInt p = pStart; p < pEnd; ++p) {
      PetscInt val;

      err = DMPlexGetLabelValue(dmMesh, name, p, &val);PYLITH_CHECK_ERROR(err);
      if (val >= 0) {
        firstPoint = p;
        break;
      } // if
    } // for
    std::string groupType = (firstPoint >= cStart && firstPoint < cEnd) ? "cell" : "vertex";
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupTypes[iGroup]), groupType);
    PetscInt numPoints;
    err = DMPlexGetStratumSize(dmMesh, name, 1, &numPoints);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(data.groupSizes[iGroup], numPoints);
    PetscIS pointIS = NULL;
    const PetscInt *points = NULL;
    const PetscInt offset = ("vertex" == groupType) ? numCells : 0;
    err = DMPlexGetStratumIS(dmMesh, name, 1, &pointIS);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    for(PetscInt p = 0; p < numPoints; ++p) {
      CPPUNIT_ASSERT_EQUAL(data.groups[index++], points[p]-offset);
    } // for
    err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
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
