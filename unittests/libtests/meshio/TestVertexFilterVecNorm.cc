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

#include "TestVertexFilterVecNorm.hh" // Implementation of class methods

#include "pylith/meshio/VertexFilterVecNorm.hh"

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestVertexFilterVecNorm );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestVertexFilterVecNorm::testConstructor(void)
{ // testConstructor
  VertexFilterVecNorm filter;
} // testConstructor

// ----------------------------------------------------------------------
// Test filter()
void
pylith::meshio::TestVertexFilterVecNorm::testFilter(void)
{ // testFilter
  const char* filename = "data/tri3.mesh";
  const int fiberDim = 2;
  const int nvertices = 4;
  const char* fieldName = "field data";
  const VectorFieldEnum fieldType = VECTOR_FIELD;
  const double fieldValues[] = {
    1.1, 1.2,
    2.1, 2.2,
    3.1, 3.2,
    4.1, 4.2
  };
  const VectorFieldEnum fieldTypeE = SCALAR_FIELD;
  const int fiberDimE = 1;
  const double fieldValuesE[] = {
    sqrt(pow(1.1, 2) + pow(1.2, 2)),
    sqrt(pow(2.1, 2) + pow(2.2, 2)),
    sqrt(pow(3.1, 2) + pow(3.2, 2)),
    sqrt(pow(4.1, 2) + pow(4.2, 2))
  };

  MeshIOAscii iohandler;
  ALE::Obj<Mesh> mesh;
  iohandler.filename(filename);
  iohandler.read(&mesh);
  CPPUNIT_ASSERT(!mesh.isNull());

  // Set vertex field
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  ALE::Obj<real_section_type> field = 
    new real_section_type(mesh->comm(), mesh->debug());
  field->setChart(mesh->getSieve()->getChart());
  field->setFiberDimension(vertices, fiberDim);
  mesh->allocate(field);

  CPPUNIT_ASSERT_EQUAL(nvertices, int(vertices->size()));

  int ipt = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt) {
    const double* values = &fieldValues[ipt*fiberDim];
    field->updatePoint(*v_iter, values);
  } // for

  VertexFilterVecNorm filter;
  VectorFieldEnum fieldTypeF = OTHER_FIELD;
  const ALE::Obj<real_section_type>& fieldF =
    filter.filter(&fieldTypeF, field, mesh);

  CPPUNIT_ASSERT_EQUAL(fieldTypeE, fieldTypeF);
  ipt = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt) {
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fieldF->getFiberDimension(*v_iter));
    const double* values = fieldF->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != values);
    const double tolerance = 1.0e-06;
    for (int i=0; i < fiberDimE; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				   values[i]/fieldValuesE[ipt*fiberDimE+i],
				   tolerance);
  } // for
} // testFilter


// End of file 
