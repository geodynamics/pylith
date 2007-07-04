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

#include "TestBoundary.hh" // Implementation of class methods

#include <Selection.hh> // USES submesh algorithms

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/utils/array.hh" // USES int_array, double_array
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "data/BoundaryDataTri3.hh" // USES BoundaryDataTri3

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestBoundary );

// ----------------------------------------------------------------------
// Test createBoundary() with 2-D triangular element.
void
pylith::faults::TestBoundary::testCreateBoundaryTri3(void)
{ // testCreateBoundaryTri3
  BoundaryDataTri3 data;
  _testCreateBoundary(data);
} // testCreateBoundaryTri3

// ----------------------------------------------------------------------
// Test createBoundary().
void
pylith::faults::TestBoundary::_testCreateBoundary(const BoundaryData& data)
{ // _testAdjustTopology
  ALE::Obj<ALE::Mesh> mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename(data.filename);
  iohandler.debug(false);
  iohandler.interpolate(false);
  iohandler.read(&mesh);

  // Extract submesh for "traction"
  const ALE::Obj<ALE::Mesh> submesh = ALE::Selection<ALE::Mesh>::submesh(mesh, mesh->getIntSection("traction"));

  CPPUNIT_ASSERT_EQUAL(data.cellDim, submesh->getDimension());

  // Check vertices
  const ALE::Obj<Mesh::label_sequence>& vertices = submesh->depthStratum(0);
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
  const ALE::Obj<sieve_type>& sieve = submesh->getSieve();
  const ALE::Obj<Mesh::label_sequence>& cells = submesh->heightStratum(1);

  const int numCells = cells->size();
  CPPUNIT_ASSERT_EQUAL(data.numCells, numCells);

  int iCell = 0;
  i = 0;
  //mesh->view(data.filename);
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
#if 0
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
#endif
} // _testCreateBoundary

// End of file 
