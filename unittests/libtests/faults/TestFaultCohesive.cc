// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFaultCohesive.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultsCohesiveKin
#include "pylith/faults/FaultCohesiveDyn.hh" // USES FaultsCohesiveDyn

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/utils/array.hh" // USES int_array, double_array
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/feassemble/Quadrature2Din3D.hh" // USES Quadrature2Din3D

#include "data/CohesiveDataLine2.hh" // USES CohesiveDataLine2
#include "data/CohesiveDataTri3.hh" // USES CohesiveDataTri3
#include "data/CohesiveDataTri3b.hh" // USES CohesiveDataTri3b
#include "data/CohesiveDataTri3c.hh" // USES CohesiveDataTri3c
#include "data/CohesiveDataTri3d.hh" // USES CohesiveDataTri3d
#include "data/CohesiveDataQuad4.hh" // USES CohesiveDataQuad4
#include "data/CohesiveDataQuad4b.hh" // USES CohesiveDataQuad4b
#include "data/CohesiveDataQuad4c.hh" // USES CohesiveDataQuad4c
#include "data/CohesiveDataQuad4d.hh" // USES CohesiveDataQuad4d
#include "data/CohesiveDataQuad4e.hh" // USES CohesiveDataQuad4e
#include "data/CohesiveDataTet4.hh" // USES CohesiveDataTet4
#include "data/CohesiveDataTet4b.hh" // USES CohesiveDataTet4b
#include "data/CohesiveDataTet4c.hh" // USES CohesiveDataTet4c
#include "data/CohesiveDataTet4d.hh" // USES CohesiveDataTet4d
#include "data/CohesiveDataHex8.hh" // USES CohesiveDataHex8
#include "data/CohesiveDataHex8b.hh" // USES CohesiveDataHex8b
#include "data/CohesiveDataHex8c.hh" // USES CohesiveDataHex8c
#include "data/CohesiveDataHex8d.hh" // USES CohesiveDataHex8d
#include "data/CohesiveDataHex8e.hh" // USES CohesiveDataHex8e
#include "data/CohesiveDataHex8f.hh" // USES CohesiveDataHex8f

#include "data/CohesiveDataLine2Lagrange.hh" // USES CohesiveDataLine2Lagrange
#include "data/CohesiveDataTri3Lagrange.hh" // USES CohesiveDataTri3Lagrange
#include "data/CohesiveDataQuad4Lagrange.hh" // USES CohesiveDataQuad4Lagrange
#include "data/CohesiveDataTet4Lagrange.hh" // USES CohesiveDataTet4Lagrange
#include "data/CohesiveDataHex8Lagrange.hh" // USES CohesiveDataHex8Lagrange

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesive );

// ----------------------------------------------------------------------
// Test adjustTopology() with 1-D line element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyLine2(void)
{ // testAdjustTopologyLine2
  CohesiveDataLine2 data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyLine2

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D triangular element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3(void)
{ // testAdjustTopologyTri3
  CohesiveDataTri3 data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyTri3

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D triangular element (edge b).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3b(void)
{ // testAdjustTopologyTri3b
  CohesiveDataTri3b data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyTri3b

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D triangular element (edge c).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3c(void)
{ // testAdjustTopologyTri3c
  CohesiveDataTri3c data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyTri3c

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D triangular element (two cohesive cells).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3d(void)
{ // testAdjustTopologyTri3d
  CohesiveDataTri3d data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyTri3d

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4(void)
{ // testAdjustTopologyQuad4
  CohesiveDataQuad4 data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyQuad4

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4b(void)
{ // testAdjustTopologyQuad4b
  CohesiveDataQuad4b data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyQuad4b

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4c(void)
{ // testAdjustTopologyQuad4c
  CohesiveDataQuad4c data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyQuad4c

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4d(void)
{ // testAdjustTopologyQuad4d
  CohesiveDataQuad4d data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyQuad4d

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4e(void)
{ // testAdjustTopologyQuad4e
  CohesiveDataQuad4e data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyQuad4e

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4(void)
{ // testAdjustTopologyTet4
  CohesiveDataTet4 data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyTet4

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element (face b).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4b(void)
{ // testAdjustTopologyTet4b
  CohesiveDataTet4b data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyTet4b

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element (face b).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4c(void)
{ // testAdjustTopologyTet4c
  CohesiveDataTet4c data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyTet4c

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element (face b).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4d(void)
{ // testAdjustTopologyTet4d
  CohesiveDataTet4d data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyTet4d

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8(void)
{ // testAdjustTopologyHex8
  CohesiveDataHex8 data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyHex8

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8b(void)
{ // testAdjustTopologyHex8b
  CohesiveDataHex8b data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyHex8b

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8c(void)
{ // testAdjustTopologyHex8c
  CohesiveDataHex8c data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyHex8c

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8d(void)
{ // testAdjustTopologyHex8d
  CohesiveDataHex8d data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyHex8d

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8e(void)
{ // testAdjustTopologyHex8e
  CohesiveDataHex8e data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyHex8e

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8f(void)
{ // testAdjustTopologyHex8f
  CohesiveDataHex8f data;
  FaultCohesiveDyn fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyHex8f

// ----------------------------------------------------------------------
// Test adjustTopology() with 1-D line element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyLine2Lagrange(void)
{ // testAdjustTopologyLine2Lagrange
  CohesiveDataLine2Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyLine2Lagrange

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D triangular element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3Lagrange(void)
{ // testAdjustTopologyTri3Lagrange
  CohesiveDataTri3Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyTri3Lagrange

// ----------------------------------------------------------------------
// Test adjustTopology() with 2-D quadrilateral element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4Lagrange(void)
{ // testAdjustTopologyQuad4Lagrange
  CohesiveDataQuad4Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyQuad4Lagrange

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D tetrahedral element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4Lagrange(void)
{ // testAdjustTopologyTet4Lagrange
  CohesiveDataTet4Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyTet4Lagrange

// ----------------------------------------------------------------------
// Test adjustTopology() with 3-D hexahedral element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8Lagrange(void)
{ // testAdjustTopologyHex8Lagrange
  CohesiveDataHex8Lagrange data;
  FaultCohesiveKin fault;
  _testAdjustTopology(&fault, data);
} // testAdjustTopologyHex8Lagrange

// ----------------------------------------------------------------------
// Test _orient1D().
void
pylith::faults::TestFaultCohesive::testOrient1D(void)
{ // testOrient1D
  double_array jacobian;
  double jacobianDet;
  double_array upDir;
  double_array orientation(1);
  
  FaultCohesive::_orient1D(&orientation, jacobian, jacobianDet, upDir);

  const int size = orientation.size();
  CPPUNIT_ASSERT_EQUAL(1, size);
  const double tolerance = 1.0e-6;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, orientation[i], tolerance);
} // testOrient1D

// ----------------------------------------------------------------------
// Test _orient2D().
void
pylith::faults::TestFaultCohesive::testOrient2D(void)
{ // testOrient2D
  const int numLocs = 2;
  const int spaceDim = 2;
  const int orientSize = 4;

  const double jacobianVals[] = {
    -1.0, 2.0,
    -0.5, 1.0
  };
  const double orientationE[] = {
    -1.0, 2.0,  -2.0, -1.0,
    -0.5, 1.0,  -1.0, -0.5
  };

  const int jacobianSize = spaceDim*(spaceDim-1);
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    double_array jacobian(&jacobianVals[iLoc*jacobianSize], jacobianSize);
    double jacobianDet;
    double_array upDir;
    double_array orientation(orientSize);

    FaultCohesive::_orient2D(&orientation, jacobian, jacobianDet, upDir);

    const int size = orientation.size();
    CPPUNIT_ASSERT_EQUAL(orientSize, size);
    const double tolerance = 1.0e-6;
    for (int i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(orientationE[iLoc*orientSize+i],
				   orientation[i], tolerance);
  } // for
} // testOrient2D

// ----------------------------------------------------------------------
// Test _orient3D().
void
pylith::faults::TestFaultCohesive::testOrient3D(void)
{ // testOrient3D
  const int numLocs = 2;
  const int spaceDim = 3;
  const int orientSize = 9;

  const double jacobianVals[] = {
    2.0,  -0.5,
    1.0,  -0.2,
    0.5,   2.0,

    -1.0,  2.0,
    -3.0, -0.2,
    -0.3,  0.3,
  };
  const double jacobianDetVals[] = {
    1.3, 0.7
  };
  const double upDirVals[] = { 0.0, 0.0, 1.0 };
  const double orientationE[] = {
    1.1654847299258313, 0.57588657243394026, 0.0, 
    -0.012145479112634533, 0.024580136299379406, 1.2997108540889502, 
    0.57575848378190342, -1.1652255028919474, 0.027417070656281111,

    0.20879249519516274, -0.66813598462452073, 0.0,
    0.65951433797689429, 0.20609823061777949, 0.11209084413599232, 
    -0.10698846644884991, -0.033433895765265592, 0.69096717914882233
  };

  double_array upDir(upDirVals, 3);
  const int jacobianSize = spaceDim*(spaceDim-1);
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    double_array jacobian(&jacobianVals[iLoc*jacobianSize], jacobianSize);
    double jacobianDet = jacobianDetVals[iLoc];
    double_array orientation(orientSize);

    FaultCohesive::_orient3D(&orientation, jacobian, jacobianDet, upDir);

    const int size = orientation.size();
    CPPUNIT_ASSERT_EQUAL(orientSize, size);
    const double tolerance = 1.0e-6;
    for (int i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(orientationE[iLoc*orientSize+i],
				   orientation[i], tolerance);
  } // for
} // testOrient3D

// ----------------------------------------------------------------------
// Test adjustTopology().
void
pylith::faults::TestFaultCohesive::_testAdjustTopology(Fault* fault,
						      const CohesiveData& data)
{ // _testAdjustTopology
  ALE::Obj<ALE::Mesh> mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename(data.filename);
  iohandler.debug(true);
  iohandler.interpolate(false);
  iohandler.read(&mesh);

  CPPUNIT_ASSERT(0 != fault);
  fault->id(1);
  fault->label("fault");
  fault->adjustTopology(mesh);

  CPPUNIT_ASSERT_EQUAL(data.cellDim, mesh->getDimension());

  // Check vertices
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const ALE::Obj<Mesh::real_section_type>& coordsField =
    mesh->getRealSection("coordinates");
  const int numVertices = vertices->size();
  CPPUNIT_ASSERT_EQUAL(data.numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(data.spaceDim, 
		       coordsField->getFiberDimension(*vertices->begin()));
  int i = 0;
  const int spaceDim = data.spaceDim;
  for(Mesh::label_sequence::iterator v_iter = 
	vertices->begin();
      v_iter != vertices->end();
      ++v_iter) {
    const Mesh::real_section_type::value_type *vertexCoords = 
      coordsField->restrictPoint(*v_iter);
    const double tolerance = 1.0e-06;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      if (data.vertices[i] < 1.0)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(data.vertices[i++], vertexCoords[iDim],
				   tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vertexCoords[iDim]/data.vertices[i++],
				   tolerance);
  } // for

  // check cells
  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

  const int numCells = cells->size();
  CPPUNIT_ASSERT_EQUAL(data.numCells, numCells);

  int iCell = 0;
  i = 0;
  mesh->view(data.filename);
  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    const int numCorners = sieve->nCone(*c_iter, mesh->depth())->size();
    CPPUNIT_ASSERT_EQUAL(data.numCorners[iCell++], numCorners);
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = 
      sieve->cone(*c_iter);
    for(sieve_type::traits::coneSequence::iterator v_iter = cone->begin();
	v_iter != cone->end();
	++v_iter)
      CPPUNIT_ASSERT_EQUAL(data.cells[i++], *v_iter);
  } // for

  // check materials
  const ALE::Obj<Mesh::label_type>& labelMaterials = 
    mesh->getLabel("material-id");
  const int idDefault = -999;
  const int size = numCells;
  int_array materialIds(size);
  i = 0;
  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter)
    materialIds[i++] = mesh->getValue(labelMaterials, *c_iter, idDefault);
  
  for (int iCell=0; iCell < numCells; ++iCell)
    CPPUNIT_ASSERT_EQUAL(data.materialIds[iCell], materialIds[iCell]);

  // Check groups
  const ALE::Obj<std::set<std::string> >& groupNames = 
    mesh->getIntSections();
  int iGroup = 0;
  int index = 0;
  for (std::set<std::string>::const_iterator name=groupNames->begin();
       name != groupNames->end();
       ++name, ++iGroup) {
    const ALE::Obj<int_section_type>& groupField = mesh->getIntSection(*name);
    CPPUNIT_ASSERT(!groupField.isNull());
    const int_section_type::chart_type& chart = groupField->getChart();
    const Mesh::point_type firstPoint = *chart.begin();
    std::string groupType = 
      (mesh->height(firstPoint) == 0) ? "cell" : "vertex";
    const int numPoints = chart.size();
    int_array points(numPoints);
    int i = 0;
    for(int_section_type::chart_type::iterator c_iter = chart.begin();
	c_iter != chart.end();
	++c_iter)
      points[i++] = *c_iter;

    CPPUNIT_ASSERT_EQUAL(std::string(data.groupNames[iGroup]), *name);
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupTypes[iGroup]), groupType);
    CPPUNIT_ASSERT_EQUAL(data.groupSizes[iGroup], numPoints);
    for (int i=0; i < numPoints; ++i)
      CPPUNIT_ASSERT_EQUAL(data.groups[index++], points[i]);
  } // for
} // _testAdjustTopology

// End of file 
