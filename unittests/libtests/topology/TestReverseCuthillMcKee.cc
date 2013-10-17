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

#include "TestReverseCuthillMcKee.hh" // Implementation of class methods

#include "pylith/topology/ReverseCuthillMcKee.hh" // USES ReverseCuthillMcKee

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

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
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestReverseCuthillMcKee );

// ----------------------------------------------------------------------
// Test reorder() with tri3 cells and no fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderTri3(void)
{ // testReorderTri3
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveTri3Level1 data;
  _testReorder(data);

  PYLITH_METHOD_END;
} // testReorderTri3

// ----------------------------------------------------------------------
// Test reorder() with level 2, tri3 cells, and one fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderTri3Fault(void)
{ // testReorderTri3Fault
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveTri3Level1Fault1 data;
  _testReorder(data);

  PYLITH_METHOD_END;
} // testReorderTri3Fault

// ----------------------------------------------------------------------
// Test reorder() with level 2, quad4 cells, and no fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderQuad4(void)
{ // testReorderQuad4
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveQuad4Level1 data;
  _testReorder(data);

  PYLITH_METHOD_END;
} // testReorderQuad4

// ----------------------------------------------------------------------
// Test reorder() with level 2, quad4 cells, and one fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderQuad4Fault(void)
{ // testReorderQuad4Fault
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveQuad4Level1Fault1 data;
  _testReorder(data);

  PYLITH_METHOD_END;
} // testReorderQuad4Fault

// ----------------------------------------------------------------------
// Test reorder() with level 2, tet4 cells, and no fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderTet4(void)
{ // testReorderTet4
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveTet4Level1 data;
  _testReorder(data);

  PYLITH_METHOD_END;
} // testReorderTet4

// ----------------------------------------------------------------------
// Test reorder() with level 2, tet4 cells, and one fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderTet4Fault(void)
{ // testReorderTet4Fault
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveTet4Level1Fault1 data;
  _testReorder(data);

  PYLITH_METHOD_END;
} // testReorderTet4Fault

// ----------------------------------------------------------------------
// Test reorder() with hex8 cells and no fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderHex8(void)
{ // testReorderHex8
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveHex8Level1 data;
  _testReorder(data);

  PYLITH_METHOD_END;
} // testReorderHex8

// ----------------------------------------------------------------------
// Test reorder() with hex8 cells and one fault.
void
pylith::topology::TestReverseCuthillMcKee::testReorderHex8Fault(void)
{ // testReorderHex8Fault
  PYLITH_METHOD_BEGIN;

  MeshDataCohesiveHex8Level1Fault1 data;
  _testReorder(data);

  PYLITH_METHOD_END;
} // testReorderHex8Fault

// ----------------------------------------------------------------------
void
pylith::topology::TestReverseCuthillMcKee::_setupMesh(Mesh* const mesh,
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
// Test reorder().
void
pylith::topology::TestReverseCuthillMcKee::_testReorder(const MeshDataCohesive& data)
{ // _testReorder
  PYLITH_METHOD_BEGIN;

  Mesh mesh(data.cellDim);
  _setupMesh(&mesh, data);

  ReverseCuthillMcKee::reorder(&mesh);
  
  // Check mesh dimension
  CPPUNIT_ASSERT_EQUAL(data.cellDim, mesh.dimension());

  const PetscDM& dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);

  // Check vertices (size only)
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  CPPUNIT_ASSERT_EQUAL(data.numVertices, verticesStratum.size());

  // Check cells (size only)
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  const PetscInt numCells = cellsStratum.size();
  CPPUNIT_ASSERT_EQUAL(data.numCells+data.numCellsCohesive, numCells);

  // Check groups
  PetscInt numGroups, pStart, pEnd;
  PetscErrorCode err;
  err = DMPlexGetChart(dmMesh, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetNumLabels(dmMesh, &numGroups);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(data.numGroups, numGroups-2); // Omit depth and material labels
  PetscInt index  = 0;
  for(PetscInt iGroup = 0; iGroup < data.numGroups; ++iGroup) {
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
  } // for

  PYLITH_METHOD_END;
} // _testReorder


// End of file 
