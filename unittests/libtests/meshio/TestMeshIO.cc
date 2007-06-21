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

#include "TestMeshIO.hh" // Implementation of class methods

#include "pylith/meshio/MeshIO.hh" // USES MeshIO
#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/utils/array.hh" // USES int_array

#include "data/MeshData.hh"

#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::logic_error

// ----------------------------------------------------------------------
// Get simple mesh for testing I/O.
ALE::Obj<ALE::Mesh>*
pylith::meshio::TestMeshIO::_createMesh(const MeshData& data)
{ // _createMesh
  // buildTopology() requires zero based index
  CPPUNIT_ASSERT(true == data.useIndexZero);

  CPPUNIT_ASSERT(0 != data.vertices);
  CPPUNIT_ASSERT(0 != data.cells);
  CPPUNIT_ASSERT(0 != data.materialIds);
  if (data.numGroups > 0) {
    CPPUNIT_ASSERT(0 != data.groups);
    CPPUNIT_ASSERT(0 != data.groupSizes);
    CPPUNIT_ASSERT(0 != data.groupNames);
    CPPUNIT_ASSERT(0 != data.groupTypes);
  } // if

  ALE::Obj<Mesh>* mesh = new ALE::Obj<Mesh>;
  CPPUNIT_ASSERT(0 != mesh);
  *mesh = new Mesh(PETSC_COMM_WORLD, data.cellDim);
  CPPUNIT_ASSERT(!mesh->isNull());
  ALE::Obj<sieve_type> sieve = new sieve_type((*mesh)->comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  // Cells and vertices
  const bool interpolate = false;
  ALE::SieveBuilder<Mesh>::buildTopology(sieve, data.cellDim, data.numCells,
			      const_cast<int*>(data.cells), data.numVertices,
					 interpolate, data.numCorners);
  (*mesh)->setSieve(sieve);
  (*mesh)->stratify();
  ALE::SieveBuilder<Mesh>::buildCoordinates(*mesh, data.spaceDim, 
					    data.vertices);

  // Material ids
  const ALE::Obj<Mesh::label_sequence>& cells = (*mesh)->heightStratum(0);
  CPPUNIT_ASSERT(!cells.isNull());
  const ALE::Obj<Mesh::label_type>& labelMaterials = 
    (*mesh)->createLabel("material-id");
  CPPUNIT_ASSERT(!labelMaterials.isNull());
  int i = 0;
  for(Mesh::label_sequence::iterator e_iter=cells->begin(); 
      e_iter != cells->end();
      ++e_iter)
    (*mesh)->setValue(labelMaterials, *e_iter, data.materialIds[i++]);

  // Groups
  for (int iGroup=0, index=0; iGroup < data.numGroups; ++iGroup) {
    const ALE::Obj<int_section_type>& groupField = 
      (*mesh)->getIntSection(data.groupNames[iGroup]);
    CPPUNIT_ASSERT(!groupField.isNull());

    MeshIO::GroupPtType type;
    const int numPoints = data.groupSizes[iGroup];
    if (0 == strcasecmp("cell", data.groupTypes[iGroup])) {
      type = MeshIO::CELL;
      for(int i=0; i < numPoints; ++i)
	groupField->setFiberDimension(data.groups[index++], 1);
    } else if (0 == strcasecmp("vertex", data.groupTypes[iGroup])) {
      type = MeshIO::VERTEX;
      const int numCells = (*mesh)->heightStratum(0)->size();
      for(int i=0; i < numPoints; ++i)
	groupField->setFiberDimension(data.groups[index++]+numCells, 1);
    } else
      throw std::logic_error("Could not parse group type.");
    (*mesh)->allocate(groupField);
  } // for
  (*mesh)->getFactory()->clear();
 
  return mesh;
} // _createMesh

// ----------------------------------------------------------------------
// Check values in mesh against data.
void
pylith::meshio::TestMeshIO::_checkVals(const ALE::Obj<Mesh>& mesh,
				       const MeshData& data)
{ // _checkVals
  // Check mesh dimension
  CPPUNIT_ASSERT_EQUAL(data.cellDim, mesh->getDimension());

  // Check vertices
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const ALE::Obj<Mesh::real_section_type>& coordsField =
    mesh->getRealSection("coordinates");
  const int numVertices = vertices->size();
  CPPUNIT_ASSERT(!vertices.isNull());
  CPPUNIT_ASSERT(!coordsField.isNull());
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
    CPPUNIT_ASSERT(0 != vertexCoords);
    const double tolerance = 1.0e-06;
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
  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

  const int numCells = cells->size();
  CPPUNIT_ASSERT_EQUAL(data.numCells, numCells);
  const int numCorners = sieve->nCone(*cells->begin(), 
				      mesh->depth())->size();
  CPPUNIT_ASSERT_EQUAL(data.numCorners, numCorners);

  const int offset = (data.useIndexZero) ? numCells : numCells-1;
  i = 0;
  for(Mesh::label_sequence::iterator e_iter = cells->begin();
      e_iter != cells->end();
      ++e_iter) {
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = 
      sieve->cone(*e_iter);
    for(sieve_type::traits::coneSequence::iterator c_iter = cone->begin();
	c_iter != cone->end();
	++c_iter)
      CPPUNIT_ASSERT_EQUAL(data.cells[i++], *c_iter-offset);
  } // for

  // check materials
  const ALE::Obj<Mesh::label_type>& labelMaterials = 
    mesh->getLabel("material-id");
  const int idDefault = -999;
  const int size = numCells;
  int_array materialIds(size);
  i = 0;
  for(Mesh::label_sequence::iterator e_iter = cells->begin();
      e_iter != cells->end();
      ++e_iter)
    materialIds[i++] = mesh->getValue(labelMaterials, *e_iter, idDefault);
  
  for (int iCell=0; iCell < numCells; ++iCell)
    CPPUNIT_ASSERT_EQUAL(data.materialIds[iCell], materialIds[iCell]);

  // Check groups
  const ALE::Obj<std::set<std::string> >& groupNames = 
    mesh->getIntSections();
  if (data.numGroups > 0) {
    CPPUNIT_ASSERT(!groupNames.isNull());
    CPPUNIT_ASSERT_EQUAL(data.numGroups, int(groupNames->size()));
  } // if
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
    const int offset = ("vertex" == groupType) ? numCells : 0;
    for(int_section_type::chart_type::iterator c_iter = chart.begin();
	c_iter != chart.end();
	++c_iter)
      points[i++] = *c_iter - offset;
    
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupNames[iGroup]), *name);
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupTypes[iGroup]), groupType);
    CPPUNIT_ASSERT_EQUAL(data.groupSizes[iGroup], numPoints);
    for (int i=0; i < numPoints; ++i)
      CPPUNIT_ASSERT_EQUAL(data.groups[index++], points[i]);
  } // for
} // _checkVals

// ----------------------------------------------------------------------
// Test debug()
void
pylith::meshio::TestMeshIO::_testDebug(MeshIO& iohandler)
{ // _testDebug
  bool debug = false;
  iohandler.debug(debug);
  CPPUNIT_ASSERT_EQUAL(debug, iohandler.debug());
  
  debug = true;
  iohandler.debug(debug);
  CPPUNIT_ASSERT_EQUAL(debug, iohandler.debug());  
} // _testDebug

// ----------------------------------------------------------------------
// Test interpolate()
void
pylith::meshio::TestMeshIO::_testInterpolate(MeshIO& iohandler)
{ // _testInterpolate
  bool interpolate = false;
  iohandler.interpolate(interpolate);
  CPPUNIT_ASSERT_EQUAL(interpolate, iohandler.interpolate());
  
  interpolate = true;
  iohandler.interpolate(interpolate);
  CPPUNIT_ASSERT_EQUAL(interpolate, iohandler.interpolate());  
} // _testInterpolate


// End of file 
