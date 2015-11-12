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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestMeshIO.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
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

  // Cells and vertices
  const bool interpolate = false;
  PetscDM dmMesh = NULL;
  PetscBool interpolateMesh = interpolate ? PETSC_TRUE : PETSC_FALSE;
  PetscInt  bound           = data.numCells*data.numCorners;
  PetscErrorCode err;


  int *cells = new int[bound];
  for (PetscInt coff = 0; coff < bound; ++coff) {cells[coff] = data.cells[coff];}
  for (PetscInt coff = 0; coff < bound; coff += data.numCorners) {
    err = DMPlexInvertCell(data.cellDim, data.numCorners, &cells[coff]);PYLITH_CHECK_ERROR(err);
  }
  err = DMPlexCreateFromCellList(_mesh->comm(), data.cellDim, data.numCells, data.numVertices, data.numCorners, interpolateMesh, cells, data.spaceDim, data.vertices, &dmMesh);PYLITH_CHECK_ERROR(err);
  delete [] cells;
  _mesh->dmMesh(dmMesh);

  // Material ids
  PetscInt cStart, cEnd;
  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    err = DMPlexSetLabelValue(dmMesh, "material-id", c, data.materialIds[c-cStart]);PYLITH_CHECK_ERROR(err);
  } // for

  // Groups
  for (int iGroup=0, index=0; iGroup < data.numGroups; ++iGroup) {
    err = DMPlexCreateLabel(dmMesh, data.groupNames[iGroup]);PYLITH_CHECK_ERROR(err);

    MeshIO::GroupPtType type;
    const int numPoints = data.groupSizes[iGroup];
    if (0 == strcasecmp("cell", data.groupTypes[iGroup])) {
      type = MeshIO::CELL;
      for(int i=0; i < numPoints; ++i, ++index) {
        err = DMPlexSetLabelValue(dmMesh, data.groupNames[iGroup], data.groups[index], 1);PYLITH_CHECK_ERROR(err);
      } // for
    } else if (0 == strcasecmp("vertex", data.groupTypes[iGroup])) {
      type = MeshIO::VERTEX;
      PetscInt numCells;
      err = DMPlexGetHeightStratum(dmMesh, 0, NULL, &numCells);PYLITH_CHECK_ERROR(err);
      for(int i=0; i < numPoints; ++i, ++index) {
        err = DMPlexSetLabelValue(dmMesh, data.groupNames[iGroup], data.groups[index]+numCells, 1);PYLITH_CHECK_ERROR(err);
      } // for
    } else
      throw std::logic_error("Could not parse group type.");
  } // for
 
    // Set coordinate system
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(data.spaceDim);
  cs.initialize();
  _mesh->coordsys(&cs);

  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(10.0);
  topology::MeshOps::nondimensionalize(_mesh, normalizer);

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
  PetscErrorCode err = 0;
  for(PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    PetscInt *closure = NULL;
    PetscInt  closureSize, numCorners = 0;

    err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    for(PetscInt p = 0; p < closureSize*2; p += 2) {
      const PetscInt point = closure[p];
      if ((point >= vStart) && (point < vEnd)) {
        closure[numCorners++] = point;
      } // if
    } // for
    err = DMPlexInvertCell(data.cellDim, numCorners, closure);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(data.numCorners, numCorners);
    for(PetscInt p = 0; p < numCorners; ++p, ++index) {
      CPPUNIT_ASSERT_EQUAL(data.cells[index], closure[p]-offset);
    } // for
    err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
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
    PetscInt numPoints, numVertices = 0;
    err = DMPlexGetStratumSize(dmMesh, name, 1, &numPoints);PYLITH_CHECK_ERROR(err);
    PetscIS pointIS = NULL;
    const PetscInt *points = NULL;
    const PetscInt offset = ("vertex" == groupType) ? numCells : 0;
    err = DMPlexGetStratumIS(dmMesh, name, 1, &pointIS);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    for(PetscInt p = 0; p < numPoints; ++p) {
      const PetscInt pStart = ("vertex" == groupType) ? vStart : cStart;
      const PetscInt pEnd   = ("vertex" == groupType) ? vEnd   : cEnd;
      if ((points[p] >= pStart) && (points[p] < pEnd)) {
        CPPUNIT_ASSERT_EQUAL(data.groups[index++], points[p]-offset);
        ++numVertices;
      }
    } // for
    CPPUNIT_ASSERT_EQUAL(data.groupSizes[iGroup], numVertices);
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
