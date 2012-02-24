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

#include "TestRefineUniform.hh" // Implementation of class methods

#include "pylith/topology/RefineUniform.hh" // USES RefineUniform

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "pylith/utils/array.hh" // USES int_array

#include "data/MeshDataCohesiveTri3Level2.hh"
#include "data/MeshDataCohesiveTri3Level2Fault1.hh"
#include "data/MeshDataCohesiveQuad4Level2.hh"
#include "data/MeshDataCohesiveQuad4Level2Fault1.hh"
#include "data/MeshDataCohesiveTet4Level2.hh"
#include "data/MeshDataCohesiveTet4Level2Fault1.hh"
#include "data/MeshDataCohesiveHex8Level2.hh"
#include "data/MeshDataCohesiveHex8Level2Fault1.hh"

#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::logic_error

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestRefineUniform );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestRefineUniform::testConstructor(void)
{ // testConstructor
  RefineUniform refiner;
} // testConstructor

// ----------------------------------------------------------------------
// Test refine() with level 2, tri3 cells, and no fault.
void
pylith::topology::TestRefineUniform::testRefineTri3Level2(void)
{ // testRefineTri3Level2
  MeshDataCohesiveTri3Level2 data;
  _testRefine(data);
} // testRefineTri3Level2

// ----------------------------------------------------------------------
// Test refine() with level 2, tri3 cells, and one fault.
void
pylith::topology::TestRefineUniform::testRefineTri3Level2Fault1(void)
{ // testRefineTri3Level2Fault1
  MeshDataCohesiveTri3Level2Fault1 data;
  _testRefine(data);
} // testRefineTri3Level2Fault1

// ----------------------------------------------------------------------
// Test refine() with level 2, quad4 cells, and no fault.
void
pylith::topology::TestRefineUniform::testRefineQuad4Level2(void)
{ // testRefineQuad4Level2
  MeshDataCohesiveQuad4Level2 data;
  _testRefine(data);
} // testRefineQuad4Level2

// ----------------------------------------------------------------------
// Test refine() with level 2, quad4 cells, and one fault.
void
pylith::topology::TestRefineUniform::testRefineQuad4Level2Fault1(void)
{ // testRefineQuad4Level2Fault1
  MeshDataCohesiveQuad4Level2Fault1 data;
  _testRefine(data);
} // testRefineQuad4Level2Fault1

// ----------------------------------------------------------------------
// Test refine() with level 2, tet4 cells, and no fault.
void
pylith::topology::TestRefineUniform::testRefineTet4Level2(void)
{ // testRefineTet4Level2
  MeshDataCohesiveTet4Level2 data;
  _testRefine(data);
} // testRefineTet4Level2

// ----------------------------------------------------------------------
// Test refine() with level 2, tet4 cells, and one fault.
void
pylith::topology::TestRefineUniform::testRefineTet4Level2Fault1(void)
{ // testRefineTet4Level2Fault1
  MeshDataCohesiveTet4Level2Fault1 data;
  _testRefine(data);
} // testRefineTet4Level2Fault1

// ----------------------------------------------------------------------
// Test refine() with level 2, hex8 cells, and no fault.
void
pylith::topology::TestRefineUniform::testRefineHex8Level2(void)
{ // testRefineHex8Level2
  MeshDataCohesiveHex8Level2 data;
  _testRefine(data);
} // testRefineHex8Level2

// ----------------------------------------------------------------------
// Test refine() with level 2, hex8 cells, and one fault.
void
pylith::topology::TestRefineUniform::testRefineHex8Level2Fault1(void)
{ // testRefineHex8Level2Fault1
  MeshDataCohesiveHex8Level2Fault1 data;
  _testRefine(data);
} // testRefineHex8Level2Fault1

// ----------------------------------------------------------------------
void
pylith::topology::TestRefineUniform::_setupMesh(Mesh* const mesh,
						const MeshDataCohesive& data)
{ // _setupMesh
  assert(0 != mesh);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(data.filename);

  iohandler.read(mesh);

  // Adjust topology if necessary.
  if (0 != data.faultA || 0 != data.faultB) {
    int firstLagrangeVertex = 0;
    int firstFaultCell = 0;

    faults::FaultCohesiveKin faultA;
    faultA.id(100);
    if (0 != data.faultA) {
      faultA.label(data.faultA);
      const int nvertices = faultA.numVerticesNoMesh(*mesh);
      firstLagrangeVertex += nvertices;
      firstFaultCell += 2*nvertices; // shadow + Lagrange vertices
    } // if

    faults::FaultCohesiveKin faultB;
    faultB.id(101);
    if (0 != data.faultB) {
      faultA.label(data.faultB);
      const int nvertices = faultB.numVerticesNoMesh(*mesh);
      firstLagrangeVertex += nvertices;
      firstFaultCell += 2*nvertices; // shadow + Lagrange vertices
    } // if
    
    int firstFaultVertex = 0;
    if (0 != data.faultA)
      faultA.adjustTopology(mesh, &firstFaultVertex, 
			    &firstLagrangeVertex, &firstFaultCell);
    if (0 != data.faultB)
      faultB.adjustTopology(mesh, &firstFaultVertex, 
			    &firstLagrangeVertex, &firstFaultCell);
  } // if
} // _setupMesh

// ----------------------------------------------------------------------
// Test refine().
void
pylith::topology::TestRefineUniform::_testRefine(const MeshDataCohesive& data)
{ // _testRefine
  typedef SieveMesh::int_section_type::chart_type chart_type;

  Mesh mesh(data.cellDim);
  _setupMesh(&mesh, data);

  RefineUniform refiner;
  Mesh newMesh(data.cellDim);
  refiner.refine(&newMesh, mesh, data.refineLevel);

  // Check mesh dimension
  CPPUNIT_ASSERT_EQUAL(data.cellDim, newMesh.dimension());

  const ALE::Obj<SieveMesh>& sieveMesh = newMesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  // Check vertices
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  const ALE::Obj<RealSection>& coordinates =
    sieveMesh->getRealSection("coordinates");
  const int numVertices = vertices->size();
  CPPUNIT_ASSERT(!vertices.isNull());
  CPPUNIT_ASSERT(!coordinates.isNull());
  CPPUNIT_ASSERT_EQUAL(data.numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(data.spaceDim, 
		       coordinates->getFiberDimension(*vertices->begin()));
  int i = 0;
  const int spaceDim = data.spaceDim;
  for(SieveMesh::label_sequence::iterator v_iter = 
	vertices->begin();
      v_iter != vertices->end();
      ++v_iter) {
    const PylithScalar* vertexCoords = coordinates->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vertexCoords);
    const PylithScalar tolerance = 1.0e-06;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      if (data.vertices[i] < 1.0) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(data.vertices[i++], vertexCoords[iDim],
				     tolerance);
      } else {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vertexCoords[iDim]/data.vertices[i++],
				     tolerance);
      }
  } // for

  // check cells
  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->heightStratum(0);

  const int numCells = cells->size();
  CPPUNIT_ASSERT_EQUAL(data.numCells+data.numCellsCohesive, numCells);

  // Normal cells
  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> pV(sieve->getMaxConeSize());
  SieveMesh::label_sequence::iterator c_iter = cells->begin();
  for(int iCell=0, i=0; iCell < data.numCells; ++iCell, ++c_iter) {
    pV.clear();
    sieve->cone(*c_iter, pV);
    const SieveMesh::point_type *cone = pV.getPoints();
    const int coneSize = pV.getSize();
    CPPUNIT_ASSERT_EQUAL(data.numCorners, coneSize);
    for(int p = 0; p < coneSize; ++p, ++i)
      CPPUNIT_ASSERT_EQUAL(data.cells[i], cone[p]);
  } // for
  // Cohesive cells
  for (int iCell=0, i=0; iCell < data.numCellsCohesive; ++iCell, ++c_iter) {
    pV.clear();
    sieve->cone(*c_iter, pV);
    const SieveMesh::point_type *cone = pV.getPoints();
    const int coneSize = pV.getSize();
    CPPUNIT_ASSERT_EQUAL(data.numCornersCohesive, coneSize);
    for(int p = 0; p < coneSize; ++p, ++i)
      CPPUNIT_ASSERT_EQUAL(data.cellsCohesive[i], cone[p]);
  } // for

  // check materials
  const ALE::Obj<SieveMesh::label_type>& labelMaterials = 
    sieveMesh->getLabel("material-id");
  const int idDefault = -999;
  const int size = numCells;
  i = 0;
  for(SieveMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    const int id = sieveMesh->getValue(labelMaterials, *c_iter, idDefault);
    CPPUNIT_ASSERT_EQUAL(data.materialIds[i++], id);
  } // for

  // Check groups
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
    for(chart_type::const_iterator c_iter = chart.begin();
	c_iter != chart.end();
	++c_iter)
      if (groupField->getFiberDimension(*c_iter))
	points[i++] = *c_iter;
    
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupNames[iGroup]), *name);
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupTypes[iGroup]), groupType);
    CPPUNIT_ASSERT_EQUAL(data.groupSizes[iGroup], numPoints);
    for (int i=0; i < numPoints; ++i)
      CPPUNIT_ASSERT_EQUAL(data.groups[index++], points[i]);
  } // for
} // _testRefine


// End of file 
