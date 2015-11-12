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

#include "TestRefineUniform.hh" // Implementation of class methods

#include "pylith/topology/RefineUniform.hh" // USES RefineUniform

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "pylith/utils/array.hh" // USES int_array

#include "data/MeshDataCohesiveTri3Level1.hh"
#include "data/MeshDataCohesiveTri3Level1Fault1.hh"
#include "data/MeshDataCohesiveQuad4Level1.hh"
#include "data/MeshDataCohesiveQuad4Level1Fault1.hh"
#include "data/MeshDataCohesiveTet4Level1.hh"
#include "data/MeshDataCohesiveTet4Level1Fault1.hh"
#include "data/MeshDataCohesiveHex8Level1.hh"
#include "data/MeshDataCohesiveHex8Level1Fault1.hh"

#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::logic_error

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestRefineUniform );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestRefineUniform::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  RefineUniform refiner;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test refine() with level 1, tri3 cells, and no fault.
void
pylith::topology::TestRefineUniform::testRefineTri3Level1(void)
{ // testRefineTri3Level1
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveTri3Level1 data;
  _testRefine(data, true);

  PYLITH_METHOD_END;
} // testRefineTri3Level1

// ----------------------------------------------------------------------
// Test refine() with level 1, tri3 cells, and one fault.
void
pylith::topology::TestRefineUniform::testRefineTri3Level1Fault1(void)
{ // testRefineTri3Level1Fault1
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveTri3Level1Fault1 data;
  _testRefine(data, true);

  PYLITH_METHOD_END;
} // testRefineTri3Level1Fault1

// ----------------------------------------------------------------------
// Test refine() with level 1, quad4 cells, and no fault.
void
pylith::topology::TestRefineUniform::testRefineQuad4Level1(void)
{ // testRefineQuad4Level1
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveQuad4Level1 data;
  _testRefine(data, false);

  PYLITH_METHOD_END;
} // testRefineQuad4Level1

// ----------------------------------------------------------------------
// Test refine() with level 1, quad4 cells, and one fault.
void
pylith::topology::TestRefineUniform::testRefineQuad4Level1Fault1(void)
{ // testRefineQuad4Level1Fault1
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveQuad4Level1Fault1 data;
  _testRefine(data, false);

  PYLITH_METHOD_END;
} // testRefineQuad4Level1Fault1

// ----------------------------------------------------------------------
// Test refine() with level 1, tet4 cells, and no fault.
void
pylith::topology::TestRefineUniform::testRefineTet4Level1(void)
{ // testRefineTet4Level1
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveTet4Level1 data;
  _testRefine(data, true);

  PYLITH_METHOD_END;
} // testRefineTet4Level1

// ----------------------------------------------------------------------
// Test refine() with level 1, tet4 cells, and one fault.
void
pylith::topology::TestRefineUniform::testRefineTet4Level1Fault1(void)
{ // testRefineTet4Level1Fault1
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveTet4Level1Fault1 data;
  _testRefine(data, true);

  PYLITH_METHOD_END;
} // testRefineTet4Level1Fault1

// ----------------------------------------------------------------------
// Test refine() with level 1, hex8 cells, and no fault.
void
pylith::topology::TestRefineUniform::testRefineHex8Level1(void)
{ // testRefineHex8Level1
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveHex8Level1 data;
  _testRefine(data, false);

  PYLITH_METHOD_END;
} // testRefineHex8Level1

// ----------------------------------------------------------------------
// Test refine() with level 1, hex8 cells, and one fault.
void
pylith::topology::TestRefineUniform::testRefineHex8Level1Fault1(void)
{ // testRefineHex8Level1Fault1
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveHex8Level1Fault1 data;
  _testRefine(data, false);

  PYLITH_METHOD_END;
} // testRefineHex8Level1Fault1

// ----------------------------------------------------------------------
void
pylith::topology::TestRefineUniform::_setupMesh(Mesh* const mesh,
						const MeshDataCohesive& data)
{ // _setupMesh
  PYLITH_METHOD_BEGIN;

  assert(mesh);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(data.filename);
  iohandler.interpolate(true);

  iohandler.read(mesh);

  // Adjust topology if necessary.
  if (data.faultA || data.faultB) {
    int firstLagrangeVertex = 0;
    int firstFaultCell = 0;

    faults::FaultCohesiveKin faultA;
    faultA.id(100);
    if (data.faultA) {
      faultA.label(data.faultA);
      const int nvertices = faultA.numVerticesNoMesh(*mesh);
      firstLagrangeVertex += nvertices;
      firstFaultCell += 2*nvertices; // shadow + Lagrange vertices
    } // if

    faults::FaultCohesiveKin faultB;
    faultB.id(101);
    if (data.faultB) {
      faultA.label(data.faultB);
      const int nvertices = faultB.numVerticesNoMesh(*mesh);
      firstLagrangeVertex += nvertices;
      firstFaultCell += 2*nvertices; // shadow + Lagrange vertices
    } // if
    
    int firstFaultVertex = 0;
    if (data.faultA) {
      faultA.adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell);
    } // if
    if (data.faultB) {
      faultB.adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell);
    } // if
  } // if

  PYLITH_METHOD_END;
} // _setupMesh

// ----------------------------------------------------------------------
// Test refine().
void
pylith::topology::TestRefineUniform::_testRefine(const MeshDataCohesive& data,
						 const bool isSimplexMesh)
{ // _testRefine
  PYLITH_METHOD_BEGIN;

  Mesh mesh(data.cellDim);
  _setupMesh(&mesh, data);

  RefineUniform refiner;
  Mesh newMesh(data.cellDim);
  refiner.refine(&newMesh, mesh, data.refineLevel);

  // Check mesh dimension
  CPPUNIT_ASSERT_EQUAL(data.cellDim, newMesh.dimension());

  const PetscDM& dmMesh = newMesh.dmMesh();CPPUNIT_ASSERT(dmMesh);

  // Check vertices
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  CPPUNIT_ASSERT_EQUAL(data.numVertices, verticesStratum.size());
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  topology::CoordsVisitor coordsVisitor(dmMesh);
  const int spaceDim = data.spaceDim;
  for (PetscInt v = vStart; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(spaceDim, coordsVisitor.sectionDof(v));
  } // for

  // Check cells
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  const PetscInt numCells = cellsStratum.size();

  CPPUNIT_ASSERT_EQUAL(data.numCells+data.numCellsCohesive, numCells);
  PetscErrorCode err;
  // Normal cells
  for(PetscInt c = cStart, index = 0; c < data.numCells; ++c) {
    PetscInt *closure = NULL;
    PetscInt closureSize, numCorners = 0;

    err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    for(PetscInt p = 0; p < closureSize*2; p += 2) {
      const PetscInt point = closure[p];
      if ((point >= vStart) && (point < vEnd)) {
        closure[numCorners++] = point;
      } // if
    } // for
    err = DMPlexInvertCell(data.cellDim, numCorners, closure);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(data.numCorners, numCorners);
    err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
  } // for

  // Cohesive cells
  for (PetscInt c = data.numCells, index = 0; c < cEnd; ++c) {
    PetscInt *closure = NULL;
    PetscInt closureSize, numCorners = 0;

    err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    for(PetscInt p = 0; p < closureSize*2; p += 2) {
      const PetscInt point = closure[p];
      if ((point >= vStart) && (point < vEnd)) {
        closure[numCorners++] = point;
      } // if
    } // for
    err = DMPlexInvertCell(data.cellDim, numCorners, closure);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(data.numCornersCohesive, numCorners);
    err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
  } // for

  // check materials
  PetscInt matId = 0;
  PetscInt matIdSum = 0; // Use sum of material ids as simple checksum.
  for(PetscInt c = cStart; c < cEnd; ++c) {
    err = DMPlexGetLabelValue(dmMesh, "material-id", c, &matId);PYLITH_CHECK_ERROR(err);
    matIdSum += matId;
  } // for
  CPPUNIT_ASSERT_EQUAL(data.matIdSum, matIdSum);

  // Check groups
  PetscInt numGroups, pStart, pEnd;
  err = DMPlexGetChart(dmMesh, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetNumLabels(dmMesh, &numGroups);PYLITH_CHECK_ERROR(err);
  PetscInt index  = 0;
  for(PetscInt iGroup = 0; iGroup < data.numGroups; ++iGroup) {
    // Omit depth, vtk, ghost and material-id labels
    // Don't know order of labels, so do brute force linear search
    bool foundLabel = false;
    int iLabel = 0;
    const char *name = NULL;
    PetscInt firstPoint = 0;

    while (iLabel < numGroups) {
      err = DMPlexGetLabelName(dmMesh, iLabel, &name);PYLITH_CHECK_ERROR(err);
      if (0 == strcmp(data.groupNames[iGroup], name)) {
        foundLabel = true;
        break;
      } else {
        ++iLabel;
      } // if/else
    } // while
    CPPUNIT_ASSERT(foundLabel);

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
    err = DMPlexGetStratumIS(dmMesh, name, 1, &pointIS);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    if (groupType == std::string("vertex")) {
      for (PetscInt p = 0; p < numPoints; ++p) {
	CPPUNIT_ASSERT((points[p] >= 0 && points[p] < cStart) || (points[p] >= cEnd));
      } // for
    } else {
      for (PetscInt p = 0; p < numPoints; ++p) {
	CPPUNIT_ASSERT(points[p] >= cStart && points[p] < cEnd);
      } // for
    } // if/else
    err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
  } // for

  PYLITH_METHOD_END;
} // _testRefine


// End of file 
