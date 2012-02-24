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

#include "TestFaultCohesive.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultsCohesiveKin
#include "pylith/faults/FaultCohesiveTract.hh" // USES FaultsCohesiveTract

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/array.hh" // USES int_array, double_array
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include "data/CohesiveDataLine2.hh" // USES CohesiveDataLine2
#include "data/CohesiveDataTri3.hh" // USES CohesiveDataTri3
#include "data/CohesiveDataTri3b.hh" // USES CohesiveDataTri3b
#include "data/CohesiveDataTri3c.hh" // USES CohesiveDataTri3c
#include "data/CohesiveDataTri3d.hh" // USES CohesiveDataTri3d
#include "data/CohesiveDataTri3e.hh" // USES CohesiveDataTri3e
#include "data/CohesiveDataTri3f.hh" // USES CohesiveDataTri3f
#include "data/CohesiveDataQuad4.hh" // USES CohesiveDataQuad4
#include "data/CohesiveDataQuad4b.hh" // USES CohesiveDataQuad4b
#include "data/CohesiveDataQuad4c.hh" // USES CohesiveDataQuad4c
#include "data/CohesiveDataQuad4d.hh" // USES CohesiveDataQuad4d
#include "data/CohesiveDataQuad4e.hh" // USES CohesiveDataQuad4e
#include "data/CohesiveDataQuad4f.hh" // USES CohesiveDataQuad4f
#include "data/CohesiveDataQuad4g.hh" // USES CohesiveDataQuad4g
#include "data/CohesiveDataQuad4h.hh" // USES CohesiveDataQuad4h
#include "data/CohesiveDataTet4.hh" // USES CohesiveDataTet4
#include "data/CohesiveDataTet4b.hh" // USES CohesiveDataTet4b
#include "data/CohesiveDataTet4c.hh" // USES CohesiveDataTet4c
#include "data/CohesiveDataTet4d.hh" // USES CohesiveDataTet4d
#include "data/CohesiveDataTet4f.hh" // USES CohesiveDataTet4f
#include "data/CohesiveDataTet4g.hh" // USES CohesiveDataTet4g
#include "data/CohesiveDataTet4h.hh" // USES CohesiveDataTet4h
#include "data/CohesiveDataTet4i.hh" // USES CohesiveDataTet4i
#include "data/CohesiveDataTet4j.hh" // USES CohesiveDataTet4j
#include "data/CohesiveDataHex8.hh" // USES CohesiveDataHex8
#include "data/CohesiveDataHex8b.hh" // USES CohesiveDataHex8b
#include "data/CohesiveDataHex8c.hh" // USES CohesiveDataHex8c
#include "data/CohesiveDataHex8d.hh" // USES CohesiveDataHex8d
#include "data/CohesiveDataHex8e.hh" // USES CohesiveDataHex8e
#include "data/CohesiveDataHex8f.hh" // USES CohesiveDataHex8f
#include "data/CohesiveDataHex8g.hh" // USES CohesiveDataHex8g
#include "data/CohesiveDataHex8h.hh" // USES CohesiveDataHex8h
#include "data/CohesiveDataHex8i.hh" // USES CohesiveDataHex8i

#include "data/CohesiveDataLine2Lagrange.hh" // USES CohesiveDataLine2Lagrange
#include "data/CohesiveDataTri3Lagrange.hh" // USES CohesiveDataTri3Lagrange
#include "data/CohesiveDataQuad4Lagrange.hh" // USES CohesiveDataQuad4Lagrange
#include "data/CohesiveDataTet4Lagrange.hh" // USES CohesiveDataTet4Lagrange
#include "data/CohesiveDataHex8Lagrange.hh" // USES CohesiveDataHex8Lagrange

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesive );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
void
pylith::faults::TestFaultCohesive::testUseFaultMesh(void)
{ // testUseFaultMesh
  FaultCohesiveTract fault;
  CPPUNIT_ASSERT(!fault._useFaultMesh);
  
  fault.useFaultMesh(true);
  CPPUNIT_ASSERT(fault._useFaultMesh);
} // testUseFaultMesh

// ----------------------------------------------------------------------
// TEMPORARY
void
pylith::faults::TestFaultCohesive::testFaultMeshFilename(void)
{ // testFaultMeshFilename
  FaultCohesiveTract fault;
  CPPUNIT_ASSERT_EQUAL(std::string("fault.inp"), fault._faultMeshFilename);
  
  const std::string filename = "SanAndreas.inp";
  fault.faultMeshFilename(filename.c_str());
  CPPUNIT_ASSERT_EQUAL(filename, fault._faultMeshFilename);
} // testUseFaultMesh

// ----------------------------------------------------------------------
// Test adjustTopology() with 1-D line element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyLine2(void)
{ // testAdjustTopologyLine2
  CohesiveDataLine2 data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyLine2

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D triangular element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3(void)
{ // testAdjustTopologyTri3
  CohesiveDataTri3 data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyTri3

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D triangular element (edge b).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3b(void)
{ // testAdjustTopologyTri3b
  CohesiveDataTri3b data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyTri3b

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D triangular element (edge c).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3c(void)
{ // testAdjustTopologyTri3c
  CohesiveDataTri3c data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyTri3c

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D triangular element (two cohesive cells).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3d(void)
{ // testAdjustTopologyTri3d
  CohesiveDataTri3d data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyTri3d

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D triangular element (two cohesive cells b).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3e(void)
{ // testAdjustTopologyTri3e
  CohesiveDataTri3e data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyTri3e

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D triangular element (vertex on fault).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3f(void)
{ // testAdjustTopologyTri3f
  CohesiveDataTri3f data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyTri3f

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4(void)
{ // testAdjustTopologyQuad4
  CohesiveDataQuad4 data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyQuad4

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4b(void)
{ // testAdjustTopologyQuad4b
  CohesiveDataQuad4b data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyQuad4b

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4c(void)
{ // testAdjustTopologyQuad4c
  CohesiveDataQuad4c data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyQuad4c

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4d(void)
{ // testAdjustTopologyQuad4d
  CohesiveDataQuad4d data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyQuad4d

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4e(void)
{ // testAdjustTopologyQuad4e
  CohesiveDataQuad4e data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyQuad4e

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4f(void)
{ // testAdjustTopologyQuad4f
  CohesiveDataQuad4f data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyQuad4f

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4g(void)
{ // testAdjustTopologyQuad4g
  CohesiveDataQuad4g data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyQuad4g

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4h(void)
{ // testAdjustTopologyQuad4h
  CohesiveDataQuad4h data;
  FaultCohesiveTract faultA;
  FaultCohesiveTract faultB;
  _testAdjustTopology(&faultA, &faultB, data, true, true);
} // testAdjustTopologyQuad4h

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4(void)
{ // testAdjustTopologyTet4
  CohesiveDataTet4 data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyTet4

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element (face b).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4b(void)
{ // testAdjustTopologyTet4b
  CohesiveDataTet4b data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyTet4b

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element (face c).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4c(void)
{ // testAdjustTopologyTet4c
  CohesiveDataTet4c data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyTet4c

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element (face d).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4d(void)
{ // testAdjustTopologyTet4d
  CohesiveDataTet4d data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyTet4d

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element (reverse cell order).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4f(void)
{ // testAdjustTopologyTet4f
  CohesiveDataTet4f data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyTet4f

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element (face g).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4g(void)
{ // testAdjustTopologyTet4g
  CohesiveDataTet4g data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyTet4g

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element (face h).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4h(void)
{ // testAdjustTopologyTet4h
  CohesiveDataTet4h data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyTet4h

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element (2 cells b).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4i(void)
{ // testAdjustTopologyTet4i
  CohesiveDataTet4i data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyTet4i

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element (vertex/edge on fault).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4j(void)
{ // testAdjustTopologyTet4j
  CohesiveDataTet4j data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyTet4j

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8(void)
{ // testAdjustTopologyHex8
  CohesiveDataHex8 data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyHex8

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8b(void)
{ // testAdjustTopologyHex8b
  CohesiveDataHex8b data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyHex8b

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8c(void)
{ // testAdjustTopologyHex8c
  CohesiveDataHex8c data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyHex8c

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8d(void)
{ // testAdjustTopologyHex8d
  CohesiveDataHex8d data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyHex8d

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8e(void)
{ // testAdjustTopologyHex8e
  CohesiveDataHex8e data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyHex8e

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8f(void)
{ // testAdjustTopologyHex8f
  CohesiveDataHex8f data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyHex8f

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element (2 cells easy).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8g(void)
{ // testAdjustTopologyHex8g
  CohesiveDataHex8g data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyHex8g

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element (2 cells difficult).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8h(void)
{ // testAdjustTopologyHex8h
  CohesiveDataHex8h data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyHex8h

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element (vertex/edge on fault).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8i(void)
{ // testAdjustTopologyHex8i
  CohesiveDataHex8i data;
  FaultCohesiveTract fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyHex8i

// ----------------------------------------------------------------------
// Test adjustTopology() with 1-D line element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyLine2Lagrange(void)
{ // testAdjustTopologyLine2Lagrange
  CohesiveDataLine2Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyLine2Lagrange

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D triangular element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3Lagrange(void)
{ // testAdjustTopologyTri3Lagrange
  CohesiveDataTri3Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyTri3Lagrange

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4Lagrange(void)
{ // testAdjustTopologyQuad4Lagrange
  CohesiveDataQuad4Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyQuad4Lagrange

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4Lagrange(void)
{ // testAdjustTopologyTet4Lagrange
  CohesiveDataTet4Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data, false);
} // testAdjustTopologyTet4Lagrange

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8Lagrange(void)
{ // testAdjustTopologyHex8Lagrange
  CohesiveDataHex8Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data, true);
} // testAdjustTopologyHex8Lagrange

// ----------------------------------------------------------------------
// Test adjustTopology().
void
pylith::faults::TestFaultCohesive::_testAdjustTopology(Fault* fault,
						       const CohesiveData& data,
						       const bool flipFault)
{ // _testAdjustTopology
  topology::Mesh mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename(data.filename);
  iohandler.debug(false);
  iohandler.interpolate(false);
  iohandler.read(&mesh);

  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh.dimension());
  cs.initialize();
  mesh.coordsys(&cs);
  mesh.nondimensionalize(normalizer);
  
  CPPUNIT_ASSERT(0 != fault);
  int firstFaultVertex    = 0;
  int firstLagrangeVertex = mesh.sieveMesh()->getIntSection("fault")->size();
  int firstFaultCell      = mesh.sieveMesh()->getIntSection("fault")->size();
  if (dynamic_cast<FaultCohesive*>(fault)->useLagrangeConstraints()) {
    firstFaultCell += mesh.sieveMesh()->getIntSection("fault")->size();
  }
  fault->id(1);
  fault->label("fault");
  fault->adjustTopology(&mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, flipFault);
  //mesh->view(data.filename);

  CPPUNIT_ASSERT_EQUAL(data.cellDim, mesh.dimension());

  // Check vertices
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const ALE::Obj<topology::Mesh::RealSection>& coordsSection =
    sieveMesh->getRealSection("coordinates");
  CPPUNIT_ASSERT(!coordsSection.isNull());
  const int numVertices = vertices->size();
  CPPUNIT_ASSERT_EQUAL(data.numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(data.spaceDim, 
		       coordsSection->getFiberDimension(*vertices->begin()));
  int i = 0;
  const int spaceDim = data.spaceDim;
  for (SieveMesh::label_sequence::iterator v_iter = 
	vertices->begin();
      v_iter != vertices->end();
      ++v_iter) {
    const double* vertexCoords = coordsSection->restrictPoint(*v_iter);
    const double tolerance = 1.0e-06;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      if (data.vertices[i] < 1.0)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(data.vertices[i++], vertexCoords[iDim],
				   tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vertexCoords[iDim]/data.vertices[i++],
				   tolerance);
  } // for

  //mesh.view("MESH");

  // check cells
  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  CPPUNIT_ASSERT(!sieve.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = sieveMesh->heightStratum(0);
  CPPUNIT_ASSERT(!cells.isNull());

  const int numCells = cells->size();
  CPPUNIT_ASSERT_EQUAL(data.numCells, numCells);

  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> pV(sieve->getMaxConeSize());
  int iCell = 0;
  i = 0;
  for(SieveMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    const int numCorners = sieveMesh->getNumCellCorners(*c_iter);
    CPPUNIT_ASSERT_EQUAL(data.numCorners[iCell++], numCorners);
    sieve->cone(*c_iter, pV);
    const SieveMesh::point_type *cone = pV.getPoints();
    for(int p = 0; p < pV.getSize(); ++p, ++i) {
      CPPUNIT_ASSERT_EQUAL(data.cells[i], cone[p]);
    }
    pV.clear();
  } // for

  // check materials
  const ALE::Obj<SieveMesh::label_type>& labelMaterials = 
    sieveMesh->getLabel("material-id");
  CPPUNIT_ASSERT(!labelMaterials.isNull());
  const int idDefault = -999;
  const int size = numCells;
  int_array materialIds(size);
  i = 0;
  for (SieveMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter)
    materialIds[i++] = sieveMesh->getValue(labelMaterials, *c_iter, idDefault);
  
  for (int iCell=0; iCell < numCells; ++iCell)
    CPPUNIT_ASSERT_EQUAL(data.materialIds[iCell], materialIds[iCell]);

  // Check groups
  const ALE::Obj<std::set<std::string> >& groupNames = 
    sieveMesh->getIntSections();
  CPPUNIT_ASSERT(!groupNames.isNull());
  int iGroup = 0;
  int index = 0;
  for (std::set<std::string>::const_iterator name=groupNames->begin();
       name != groupNames->end();
       ++name, ++iGroup) {
    const ALE::Obj<topology::Mesh::IntSection>& groupField =
      sieveMesh->getIntSection(*name);
    CPPUNIT_ASSERT(!groupField.isNull());
    const topology::Mesh::IntSection::chart_type& chart = groupField->getChart();
    SieveMesh::point_type firstPoint;
    for (topology::Mesh::IntSection::chart_type::const_iterator c_iter = chart.begin(); c_iter != chart.end(); ++c_iter) {
      if (groupField->getFiberDimension(*c_iter)) {firstPoint = *c_iter; break;}
    }
    std::string groupType = 
      (sieveMesh->height(firstPoint) == 0) ? "cell" : "vertex";
    const int numPoints = groupField->size();
    int_array points(numPoints);
    int i = 0;
    for (topology::Mesh::IntSection::chart_type::const_iterator c_iter = chart.begin(); c_iter != chart.end(); ++c_iter) {
      if (groupField->getFiberDimension(*c_iter)) points[i++] = *c_iter;
    }

    CPPUNIT_ASSERT_EQUAL(std::string(data.groupNames[iGroup]), *name);
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupTypes[iGroup]), groupType);
    CPPUNIT_ASSERT_EQUAL(data.groupSizes[iGroup], numPoints);
    for (int i=0; i < numPoints; ++i)
      CPPUNIT_ASSERT_EQUAL(data.groups[index++], points[i]);
  } // for
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
  topology::Mesh mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename(data.filename);
  iohandler.debug(false);
  iohandler.interpolate(false);
  iohandler.read(&mesh);

  CPPUNIT_ASSERT(0 != faultA);
  CPPUNIT_ASSERT(0 != faultB);
  int firstFaultVertex    = 0;
  int firstLagrangeVertex = mesh.sieveMesh()->getIntSection("faultA")->size() + mesh.sieveMesh()->getIntSection("faultB")->size();
  int firstFaultCell      = mesh.sieveMesh()->getIntSection("faultA")->size() + mesh.sieveMesh()->getIntSection("faultB")->size();
  if (dynamic_cast<FaultCohesive*>(faultA)->useLagrangeConstraints()) {
    firstFaultCell += mesh.sieveMesh()->getIntSection("faultA")->size();
  }
  if (dynamic_cast<FaultCohesive*>(faultB)->useLagrangeConstraints()) {
    firstFaultCell += mesh.sieveMesh()->getIntSection("faultB")->size();
  }

  faultA->id(1);
  faultA->label("faultA");
  faultA->adjustTopology(&mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, flipFaultA);

  faultB->id(2);
  faultB->label("faultB");
  faultB->adjustTopology(&mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, flipFaultB);

  //sieveMesh->view(data.filename);
  CPPUNIT_ASSERT_EQUAL(data.cellDim, mesh.dimension());

  // Check vertices
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const ALE::Obj<topology::Mesh::RealSection>& coordsSection =
    sieveMesh->getRealSection("coordinates");
  CPPUNIT_ASSERT(!coordsSection.isNull());
  const int numVertices = vertices->size();
  CPPUNIT_ASSERT_EQUAL(data.numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(data.spaceDim, 
		       coordsSection->getFiberDimension(*vertices->begin()));
  int i = 0;
  const int spaceDim = data.spaceDim;
  for(SieveMesh::label_sequence::iterator v_iter = 
	vertices->begin();
      v_iter != vertices->end();
      ++v_iter) {
    const double* coordsVertex = coordsSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != coordsVertex);
    const double tolerance = 1.0e-06;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      if (data.vertices[i] < 1.0)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(data.vertices[i++], coordsVertex[iDim],
				   tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, coordsVertex[iDim]/data.vertices[i++],
				   tolerance);
  } // for

  // check cells
  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  CPPUNIT_ASSERT(!sieve.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = sieveMesh->heightStratum(0);
  CPPUNIT_ASSERT(!cells.isNull());

  const int numCells = cells->size();
  CPPUNIT_ASSERT_EQUAL(data.numCells, numCells);

  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> pV(sieve->getMaxConeSize());
  int iCell = 0;
  i = 0;
  //sieveMesh->view(data.filename);
  for(SieveMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    const int numCorners = sieveMesh->getNumCellCorners(*c_iter);
    CPPUNIT_ASSERT_EQUAL(data.numCorners[iCell++], numCorners);
    sieve->cone(*c_iter, pV);
    const SieveMesh::point_type *cone = pV.getPoints();
    for(int p = 0; p < pV.getSize(); ++p, ++i) {
      CPPUNIT_ASSERT_EQUAL(data.cells[i], cone[p]);
    }
    pV.clear();
  } // for

  // check materials
  const ALE::Obj<SieveMesh::label_type>& labelMaterials = 
    sieveMesh->getLabel("material-id");
  CPPUNIT_ASSERT(!labelMaterials.isNull());
  const int idDefault = -999;
  const int size = numCells;
  int_array materialIds(size);
  i = 0;
  for(SieveMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter)
    materialIds[i++] = sieveMesh->getValue(labelMaterials, *c_iter, idDefault);
  
  for (int iCell=0; iCell < numCells; ++iCell)
    CPPUNIT_ASSERT_EQUAL(data.materialIds[iCell], materialIds[iCell]);

  // Check groups
  const ALE::Obj<std::set<std::string> >& groupNames = 
    sieveMesh->getIntSections();
  CPPUNIT_ASSERT(!groupNames.isNull());
  int iGroup = 0;
  int index = 0;
  for (std::set<std::string>::const_iterator name=groupNames->begin();
       name != groupNames->end();
       ++name, ++iGroup) {
    const ALE::Obj<topology::Mesh::IntSection>& groupField = 
      sieveMesh->getIntSection(*name);
    CPPUNIT_ASSERT(!groupField.isNull());
    const topology::Mesh::IntSection::chart_type& chart = groupField->getChart();
    SieveMesh::point_type firstPoint;
    for(topology::Mesh::IntSection::chart_type::const_iterator c_iter = chart.begin(); c_iter != chart.end(); ++c_iter) {
      if (groupField->getFiberDimension(*c_iter)) {firstPoint = *c_iter; break;}
    }
    std::string groupType = 
      (sieveMesh->height(firstPoint) == 0) ? "cell" : "vertex";
    const int numPoints = groupField->size();
    int_array points(numPoints);
    int i = 0;
    for(topology::Mesh::IntSection::chart_type::const_iterator c_iter = chart.begin(); c_iter != chart.end(); ++c_iter) {
      if (groupField->getFiberDimension(*c_iter)) points[i++] = *c_iter;
    }

    CPPUNIT_ASSERT_EQUAL(std::string(data.groupNames[iGroup]), *name);
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupTypes[iGroup]), groupType);
    CPPUNIT_ASSERT_EQUAL(data.groupSizes[iGroup], numPoints);
    for (int i=0; i < numPoints; ++i)
      CPPUNIT_ASSERT_EQUAL(data.groups[index++], points[i]);
  } // for
} // _testAdjustTopology


// End of file 
