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

#include "TestFaultCohesive.hh" // Implementation of class methods

#include "TestFaultMesh.hh" // USES createFaultMesh()

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultsCohesiveKin
#include "pylith/faults/FaultCohesiveTract.hh" // USES FaultsCohesiveTract

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/utils/array.hh" // USES int_array, scalar_array
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesive );

// ----------------------------------------------------------------------
void
pylith::faults::TestFaultCohesive::testUseFaultMesh(void)
{ // testUseFaultMesh
  PYLITH_METHOD_BEGIN;

  FaultCohesiveTract fault;
  CPPUNIT_ASSERT(!fault._useFaultMesh);
  
  fault.useFaultMesh(true);
  CPPUNIT_ASSERT(fault._useFaultMesh);

  PYLITH_METHOD_END;
} // testUseFaultMesh

// ----------------------------------------------------------------------
#include "data/CohesiveDataLine2.hh" // USES CohesiveDataLine2

// Test adjustTopology() with 1-D line element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyLine2(void)
{ // testAdjustTopologyLine2
  PYLITH_METHOD_BEGIN;

  CohesiveDataLine2 data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyLine2

// ----------------------------------------------------------------------
#include "data/CohesiveDataTri3.hh" // USES CohesiveDataTri3

// Test adjustTopology() with 2-D triangular element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3(void)
{ // testAdjustTopologyTri3
  PYLITH_METHOD_BEGIN;

  CohesiveDataTri3 data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyTri3

// ----------------------------------------------------------------------
#include "data/CohesiveDataTri3b.hh" // USES CohesiveDataTri3b

// Test adjustTopology() with 2-D triangular element (edge b).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3b(void)
{ // testAdjustTopologyTri3b
  PYLITH_METHOD_BEGIN;

  CohesiveDataTri3b data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyTri3b

// ----------------------------------------------------------------------
#include "data/CohesiveDataTri3c.hh" // USES CohesiveDataTri3c

// Test adjustTopology() with 2-D triangular element (edge c).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3c(void)
{ // testAdjustTopologyTri3c
  PYLITH_METHOD_BEGIN;

  CohesiveDataTri3c data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyTri3c

// ----------------------------------------------------------------------
#include "data/CohesiveDataTri3d.hh" // USES CohesiveDataTri3d

// Test adjustTopology() with 2-D triangular element (two cohesive cells).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3d(void)
{ // testAdjustTopologyTri3d
  PYLITH_METHOD_BEGIN;

  CohesiveDataTri3d data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyTri3d

// ----------------------------------------------------------------------
#include "data/CohesiveDataTri3e.hh" // USES CohesiveDataTri3e

// Test adjustTopology() with 2-D triangular element (two cohesive cells b).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3e(void)
{ // testAdjustTopologyTri3e
  PYLITH_METHOD_BEGIN;

  CohesiveDataTri3e data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyTri3e

// ----------------------------------------------------------------------
#include "data/CohesiveDataTri3f.hh" // USES CohesiveDataTri3f

// Test adjustTopology() with 2-D triangular element (vertex on fault).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3f(void)
{ // testAdjustTopologyTri3f
  PYLITH_METHOD_BEGIN;

  CohesiveDataTri3f data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyTri3f

// ----------------------------------------------------------------------
#include "data/CohesiveDataQuad4.hh" // USES CohesiveDataQuad4

// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4(void)
{ // testAdjustTopologyQuad4
  PYLITH_METHOD_BEGIN;

  CohesiveDataQuad4 data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyQuad4

// ----------------------------------------------------------------------
#include "data/CohesiveDataQuad4b.hh" // USES CohesiveDataQuad4b

// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4b(void)
{ // testAdjustTopologyQuad4b
  PYLITH_METHOD_BEGIN;

  CohesiveDataQuad4b data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyQuad4b

// ----------------------------------------------------------------------
#include "data/CohesiveDataQuad4c.hh" // USES CohesiveDataQuad4c

// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4c(void)
{ // testAdjustTopologyQuad4c
  PYLITH_METHOD_BEGIN;

  CohesiveDataQuad4c data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyQuad4c

// ----------------------------------------------------------------------
#include "data/CohesiveDataQuad4d.hh" // USES CohesiveDataQuad4d

// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4d(void)
{ // testAdjustTopologyQuad4d
  PYLITH_METHOD_BEGIN;

  CohesiveDataQuad4d data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyQuad4d

// ----------------------------------------------------------------------
#include "data/CohesiveDataQuad4e.hh" // USES CohesiveDataQuad4e

// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4e(void)
{ // testAdjustTopologyQuad4e
  PYLITH_METHOD_BEGIN;

  CohesiveDataQuad4e data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyQuad4e

// ----------------------------------------------------------------------
#include "data/CohesiveDataQuad4f.hh" // USES CohesiveDataQuad4f

// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4f(void)
{ // testAdjustTopologyQuad4f
  PYLITH_METHOD_BEGIN;

  CohesiveDataQuad4f data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyQuad4f

// ----------------------------------------------------------------------
#include "data/CohesiveDataQuad4g.hh" // USES CohesiveDataQuad4g

// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4g(void)
{ // testAdjustTopologyQuad4g
  PYLITH_METHOD_BEGIN;

  CohesiveDataQuad4g data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyQuad4g

// ----------------------------------------------------------------------
#include "data/CohesiveDataQuad4h.hh" // USES CohesiveDataQuad4h

// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4h(void)
{ // testAdjustTopologyQuad4h
  PYLITH_METHOD_BEGIN;

  CohesiveDataQuad4h data;
  FaultCohesiveTract faultA;
  FaultCohesiveTract faultB;
  _testAdjustTopology(&faultA, &faultB, data, true, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyQuad4h

// ----------------------------------------------------------------------
#include "data/CohesiveDataTet4.hh" // USES CohesiveDataTet4

// Test adjustTopology() with 3-D tetrahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4(void)
{ // testAdjustTopologyTet4
  PYLITH_METHOD_BEGIN;

  CohesiveDataTet4 data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyTet4

// ----------------------------------------------------------------------
#include "data/CohesiveDataTet4b.hh" // USES CohesiveDataTet4b

// Test adjustTopology() with 3-D tetrahedral element (face b).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4b(void)
{ // testAdjustTopologyTet4b
  PYLITH_METHOD_BEGIN;

  CohesiveDataTet4b data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyTet4b

// ----------------------------------------------------------------------
#include "data/CohesiveDataTet4c.hh" // USES CohesiveDataTet4c

// Test adjustTopology() with 3-D tetrahedral element (face c).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4c(void)
{ // testAdjustTopologyTet4c
  PYLITH_METHOD_BEGIN;

  CohesiveDataTet4c data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyTet4c

// ----------------------------------------------------------------------
#include "data/CohesiveDataTet4d.hh" // USES CohesiveDataTet4d

// Test adjustTopology() with 3-D tetrahedral element (face d).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4d(void)
{ // testAdjustTopologyTet4d
  PYLITH_METHOD_BEGIN;

  CohesiveDataTet4d data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyTet4d

// ----------------------------------------------------------------------
#include "data/CohesiveDataTet4f.hh" // USES CohesiveDataTet4f

// Test adjustTopology() with 3-D tetrahedral element (reverse cell order).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4f(void)
{ // testAdjustTopologyTet4f
  PYLITH_METHOD_BEGIN;

  CohesiveDataTet4f data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyTet4f

// ----------------------------------------------------------------------
#include "data/CohesiveDataTet4g.hh" // USES CohesiveDataTet4g

// Test adjustTopology() with 3-D tetrahedral element (face g).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4g(void)
{ // testAdjustTopologyTet4g
  PYLITH_METHOD_BEGIN;

  CohesiveDataTet4g data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyTet4g

// ----------------------------------------------------------------------
#include "data/CohesiveDataTet4h.hh" // USES CohesiveDataTet4h

// Test adjustTopology() with 3-D tetrahedral element (face h).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4h(void)
{ // testAdjustTopologyTet4h
  PYLITH_METHOD_BEGIN;

  CohesiveDataTet4h data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyTet4h

// ----------------------------------------------------------------------
#include "data/CohesiveDataTet4i.hh" // USES CohesiveDataTet4i

// Test adjustTopology() with 3-D tetrahedral element (2 cells b).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4i(void)
{ // testAdjustTopologyTet4i
  PYLITH_METHOD_BEGIN;

  CohesiveDataTet4i data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyTet4i

// ----------------------------------------------------------------------
#include "data/CohesiveDataTet4j.hh" // USES CohesiveDataTet4j

// Test adjustTopology() with 3-D tetrahedral element (vertex/edge on fault).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4j(void)
{ // testAdjustTopologyTet4j
  PYLITH_METHOD_BEGIN;

  CohesiveDataTet4j data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyTet4j

// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8.hh" // USES CohesiveDataHex8

// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8(void)
{ // testAdjustTopologyHex8
  PYLITH_METHOD_BEGIN;

  CohesiveDataHex8 data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyHex8

// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8b.hh" // USES CohesiveDataHex8b

// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8b(void)
{ // testAdjustTopologyHex8b
  PYLITH_METHOD_BEGIN;

  CohesiveDataHex8b data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyHex8b

// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8c.hh" // USES CohesiveDataHex8c

// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8c(void)
{ // testAdjustTopologyHex8c
  PYLITH_METHOD_BEGIN;

  CohesiveDataHex8c data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyHex8c

// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8d.hh" // USES CohesiveDataHex8d

// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8d(void)
{ // testAdjustTopologyHex8d
  PYLITH_METHOD_BEGIN;

  CohesiveDataHex8d data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyHex8d

// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8e.hh" // USES CohesiveDataHex8e

// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8e(void)
{ // testAdjustTopologyHex8e
  PYLITH_METHOD_BEGIN;

  CohesiveDataHex8e data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyHex8e

// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8f.hh" // USES CohesiveDataHex8f

// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8f(void)
{ // testAdjustTopologyHex8f
  PYLITH_METHOD_BEGIN;

  CohesiveDataHex8f data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyHex8f

// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8g.hh" // USES CohesiveDataHex8g

// Test adjustTopology() with 3-D hexahedral element (2 cells easy).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8g(void)
{ // testAdjustTopologyHex8g
  PYLITH_METHOD_BEGIN;

  CohesiveDataHex8g data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyHex8g

// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8h.hh" // USES CohesiveDataHex8h

// Test adjustTopology() with 3-D hexahedral element (2 cells difficult).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8h(void)
{ // testAdjustTopologyHex8h
  PYLITH_METHOD_BEGIN;

  CohesiveDataHex8h data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyHex8h

// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8i.hh" // USES CohesiveDataHex8i

// Test adjustTopology() with 3-D hexahedral element (vertex/edge on fault).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8i(void)
{ // testAdjustTopologyHex8i
  PYLITH_METHOD_BEGIN;

  CohesiveDataHex8i data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyHex8i

// ----------------------------------------------------------------------
#include "data/CohesiveDataLine2Lagrange.hh" // USES CohesiveDataLine2Lagrange

// Test adjustTopology() with 1-D line element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyLine2Lagrange(void)
{ // testAdjustTopologyLine2Lagrange
  PYLITH_METHOD_BEGIN;

  CohesiveDataLine2Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyLine2Lagrange

// ----------------------------------------------------------------------
#include "data/CohesiveDataTri3Lagrange.hh" // USES CohesiveDataTri3Lagrange

// Test adjustTopology() with 2-D triangular element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3Lagrange(void)
{ // testAdjustTopologyTri3Lagrange
  PYLITH_METHOD_BEGIN;

  CohesiveDataTri3Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyTri3Lagrange

// ----------------------------------------------------------------------
#include "data/CohesiveDataQuad4Lagrange.hh" // USES CohesiveDataQuad4Lagrange

// Test adjustTopology() with 2-D quadrilateral element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4Lagrange(void)
{ // testAdjustTopologyQuad4Lagrange
  PYLITH_METHOD_BEGIN;

  CohesiveDataQuad4Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyQuad4Lagrange

// ----------------------------------------------------------------------
#include "data/CohesiveDataTet4Lagrange.hh" // USES CohesiveDataTet4Lagrange

// Test adjustTopology() with 3-D tetrahedral element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4Lagrange(void)
{ // testAdjustTopologyTet4Lagrange
  PYLITH_METHOD_BEGIN;

  CohesiveDataTet4Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data, false);

  PYLITH_METHOD_END;
} // testAdjustTopologyTet4Lagrange

// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8Lagrange.hh" // USES CohesiveDataHex8Lagrange

// Test adjustTopology() with 3-D hexahedral element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8Lagrange(void)
{ // testAdjustTopologyHex8Lagrange
  PYLITH_METHOD_BEGIN;

  CohesiveDataHex8Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data, true);

  PYLITH_METHOD_END;
} // testAdjustTopologyHex8Lagrange

// ----------------------------------------------------------------------
// Test adjustTopology().
void
pylith::faults::TestFaultCohesive::_testAdjustTopology(Fault* fault,
						       const CohesiveData& data,
						       const bool flipFault)
{ // _testAdjustTopology
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(fault);

  topology::Mesh mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename(data.filename);
  iohandler.debug(false);
  iohandler.interpolate(false);
  iohandler.read(&mesh);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh.dimension());
  cs.initialize();
  mesh.coordsys(&cs);

  spatialdata::units::Nondimensional normalizer;
  topology::MeshOps::nondimensionalize(&mesh, normalizer);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscInt firstFaultVertex = 0;
  PetscInt firstLagrangeVertex = 0, firstFaultCell = 0;
  PetscErrorCode err;

  err = DMPlexGetStratumSize(dmMesh, "fault", 1, &firstLagrangeVertex);PYLITH_CHECK_ERROR(err);
  firstFaultCell = firstLagrangeVertex;
  if (dynamic_cast<FaultCohesive*>(fault)->useLagrangeConstraints()) {
    firstFaultCell += firstLagrangeVertex;
  } // if
  fault->id(1);
  fault->label("fault");
  fault->adjustTopology(&mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, flipFault);
  //mesh->view(data.filename);

  CPPUNIT_ASSERT_EQUAL(data.cellDim, mesh.dimension());
  dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);

  // Check vertices
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  CPPUNIT_ASSERT_EQUAL(data.numVertices, verticesStratum.size());
  
  topology::CoordsVisitor coordsVisitor(dmMesh);
  const PetscScalar* coordsArray = coordsVisitor.localArray();CPPUNIT_ASSERT(coordsArray);
  const PetscInt spaceDim    = data.spaceDim;

  const PylithScalar tolerance = 1.0e-06;
  for (PetscInt v = vStart, index = 0; v < vEnd; ++v) {
    const PetscInt off = coordsVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, coordsVisitor.sectionDof(v));
    for (PetscInt d = 0; d < spaceDim; ++d, ++index)
      if (fabs(data.vertices[index]) < 1.0)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(data.vertices[index], coordsArray[off+d], tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, coordsArray[off+d]/data.vertices[index], tolerance);
  } // for

  // check cells
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  CPPUNIT_ASSERT_EQUAL(data.numCells, cellsStratum.size());
  for (PetscInt c = cStart, cell = 0, i = 0; c < cEnd; ++c, ++cell) {
    PetscInt *closure = PETSC_NULL;
    PetscInt  closureSize, numCorners = 0;

    err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    for(PetscInt p = 0; p < closureSize*2; p += 2) {
      const PetscInt point = closure[p];
      if ((point >= vStart) && (point < vEnd)) {
        closure[numCorners++] = point;
      } // if
    } // for
    err = DMPlexInvertCell(data.cellDim, numCorners, closure);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(data.numCorners[cell], numCorners);
    for (PetscInt p = 0; p < numCorners; ++p, ++i) {
      CPPUNIT_ASSERT_EQUAL(data.cells[i], closure[p]);
    } // for
    err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
  } // for

  // check materials
  PetscDMLabel labelMaterials = NULL;
  err = DMPlexGetLabel(dmMesh, "material-id", &labelMaterials);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT(labelMaterials);
  const PetscInt idDefault = -999;
  for (PetscInt c = cStart, cell = 0; c < cEnd; ++c, ++cell) {
    PetscInt value;

    err = DMLabelGetValue(labelMaterials, c, &value);PYLITH_CHECK_ERROR(err);
    if (value == -1)
      value = idDefault;
    CPPUNIT_ASSERT_EQUAL(data.materialIds[cell], value);
  }  // for

  // Check groups
  PetscInt numLabels;
  err = DMPlexGetNumLabels(dmMesh, &numLabels);PYLITH_CHECK_ERROR(err);
  for (PetscInt l = numLabels-1, i = 0, index = 0; l >= 0; --l) {
    PetscDMLabel label = NULL;
    PetscIS is = NULL;
    const PetscInt *points = NULL;
    PetscInt numPoints = 0, depth = 0;
    const char *name = NULL;
    std::string skipA = "depth";
    std::string skipB = "material-id";

    err = DMPlexGetLabelName(dmMesh, l, &name);PYLITH_CHECK_ERROR(err);
    if (std::string(name) == skipA) continue;
    if (std::string(name) == skipB) continue;
    err = DMPlexGetLabel(dmMesh, name, &label);PYLITH_CHECK_ERROR(err);CPPUNIT_ASSERT(label);
    err = DMLabelGetStratumIS(label, 1, &is);PYLITH_CHECK_ERROR(err);
    err = ISGetLocalSize(is, &numPoints);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(is, &points);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetLabelValue(dmMesh, "depth", points[0], &depth);PYLITH_CHECK_ERROR(err);
    std::string groupType = depth ? "cell" : "vertex";

    CPPUNIT_ASSERT_EQUAL(std::string(data.groupNames[i]), std::string(name));
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupTypes[i]), groupType);
    CPPUNIT_ASSERT_EQUAL(data.groupSizes[i], numPoints);
    for (PetscInt p = 0; p < numPoints; ++p, ++index)
      CPPUNIT_ASSERT_EQUAL(data.groups[index], points[p]);
    err = ISRestoreIndices(is, &points);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&is);PYLITH_CHECK_ERROR(err);
    ++i;
  } // for

  PYLITH_METHOD_END;
} // _testAdjustTopology

// ----------------------------------------------------------------------
// Test adjustTopology().
void
pylith::faults::TestFaultCohesive::_testAdjustTopology(Fault* faultA,
						       Fault* faultB,
						       const CohesiveData& data,
						       const bool flipFaultA,
						       const bool flipFaultB)
{ // _testAdjustTopology
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(faultA);
  CPPUNIT_ASSERT(faultB);

  topology::Mesh mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename(data.filename);
  iohandler.debug(false);
  iohandler.interpolate(false);
  iohandler.read(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscInt firstFaultVertex = 0;
  PetscInt sizeA = 0, sizeB = 0;
  PetscErrorCode err;

  err = DMPlexGetStratumSize(dmMesh, "faultA", 1, &sizeA);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetStratumSize(dmMesh, "faultB", 1, &sizeB);PYLITH_CHECK_ERROR(err);
  PetscInt firstLagrangeVertex = sizeA + sizeB;
  PetscInt firstFaultCell = sizeA + sizeB;
  if (dynamic_cast<FaultCohesive*>(faultA)->useLagrangeConstraints()) {
    firstFaultCell += sizeA;
  } // if
  if (dynamic_cast<FaultCohesive*>(faultB)->useLagrangeConstraints()) {
    firstFaultCell += sizeB;
  } // if

  faultA->id(1);
  faultA->label("faultA");
  faultA->adjustTopology(&mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, flipFaultA);

  faultB->id(2);
  faultB->label("faultB");
  faultB->adjustTopology(&mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, flipFaultB);

  CPPUNIT_ASSERT_EQUAL(data.cellDim, mesh.dimension());
  dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);

  // Check vertices
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  CPPUNIT_ASSERT_EQUAL(data.numVertices, verticesStratum.size());
  
  topology::CoordsVisitor coordsVisitor(dmMesh);
  const PetscScalar* coordsArray = coordsVisitor.localArray();CPPUNIT_ASSERT(coordsArray);
  const PetscInt spaceDim = data.spaceDim;

  const PylithScalar tolerance = 1.0e-06;
  for (PetscInt v = vStart, index = 0; v < vEnd; ++v) {
    const PetscInt off = coordsVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, coordsVisitor.sectionDof(v));
    for (PetscInt d = 0; d < spaceDim; ++d, ++index)
      if (fabs(data.vertices[index]) < 1.0)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(data.vertices[index], coordsArray[off+d], tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, coordsArray[off+d]/data.vertices[index], tolerance);
  } // for

  // check cells
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  CPPUNIT_ASSERT_EQUAL(data.numCells, cellsStratum.size());
  for (PetscInt c = cStart, cell = 0, i = 0; c < cEnd; ++c, ++cell) {
    const PetscInt *cone = NULL;
    PetscInt coneSize = 0;
    err = DMPlexGetConeSize(dmMesh, c, &coneSize);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetCone(dmMesh, c, &cone);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(data.numCorners[cell], coneSize);
    for (PetscInt p = 0; p < coneSize; ++p, ++i) {
      CPPUNIT_ASSERT_EQUAL(data.cells[i], cone[p]);
    } // for
  } // for

  // check materials
  PetscDMLabel labelMaterials = NULL;
  err = DMPlexGetLabel(dmMesh, "material-id", &labelMaterials);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT(labelMaterials);
  const PetscInt idDefault = -999;
  for (PetscInt c = cStart, cell = 0; c < cEnd; ++c, ++cell) {
    PetscInt value;

    err = DMLabelGetValue(labelMaterials, c, &value);PYLITH_CHECK_ERROR(err);
    if (value == -1)
      value = idDefault;
    CPPUNIT_ASSERT_EQUAL(data.materialIds[cell], value);
  }  

  // Check groups
  PetscInt numLabels;
  err = DMPlexGetNumLabels(dmMesh, &numLabels);PYLITH_CHECK_ERROR(err);
  for (PetscInt l = 0, i = 0, index = 0; l < numLabels; ++l) {
    PetscDMLabel label = NULL;
    PetscIS is = NULL;
    const PetscInt *points = NULL;
    PetscInt numPoints = 0, depth = 0;
    const char *name = NULL;
    std::string skipA = "depth";
    std::string skipB = "material-id";

    err = DMPlexGetLabelName(dmMesh, l, &name);PYLITH_CHECK_ERROR(err);
    if (std::string(name) == skipA) continue;
    if (std::string(name) == skipB) continue;
    err = DMPlexGetLabel(dmMesh, name, &label);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT(label);
    err = DMLabelGetStratumIS(label, 1, &is);PYLITH_CHECK_ERROR(err);
    err = ISGetLocalSize(is, &numPoints);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(is, &points);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetLabelValue(dmMesh, "depth", points[0], &depth);PYLITH_CHECK_ERROR(err);
    std::string groupType = depth ? "cell" : "vertex";

    CPPUNIT_ASSERT_EQUAL(std::string(data.groupNames[i]), std::string(name));
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupTypes[i]), groupType);
    CPPUNIT_ASSERT_EQUAL(data.groupSizes[i], numPoints);
    for (PetscInt p = 0; p < numPoints; ++p, ++index)
      CPPUNIT_ASSERT_EQUAL(data.groups[index], points[p]);
    err = ISRestoreIndices(is, &points);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&is);PYLITH_CHECK_ERROR(err);
    ++i;
  } // for

  PYLITH_METHOD_END;
} // _testAdjustTopology


// End of file 
