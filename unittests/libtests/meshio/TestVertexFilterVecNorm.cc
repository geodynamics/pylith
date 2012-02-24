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

#include "TestVertexFilterVecNorm.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/meshio/VertexFilterVecNorm.hh"

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/Field.hh" // USES Field

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestVertexFilterVecNorm );

// ----------------------------------------------------------------------
typedef pylith::topology::Field<pylith::topology::Mesh> MeshField;

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestVertexFilterVecNorm::testConstructor(void)
{ // testConstructor
  VertexFilterVecNorm<MeshField> filter;
} // testConstructor

// ----------------------------------------------------------------------
// Test filter()
void
pylith::meshio::TestVertexFilterVecNorm::testFilter(void)
{ // testFilter
  typedef pylith::topology::Mesh::SieveMesh SieveMesh;
  typedef pylith::topology::Mesh::RealSection RealSection;

  const char* filename = "data/tri3.mesh";
  const int fiberDim = 2;
  const int nvertices = 4;
  const std::string label = "field data";
  const topology::FieldBase::VectorFieldEnum fieldType = 
    topology::FieldBase::VECTOR;
  const PylithScalar fieldValues[] = {
    1.1, 1.2,
    2.1, 2.2,
    3.1, 3.2,
    4.1, 4.2
  };
  const topology::FieldBase::VectorFieldEnum fieldTypeE = 
    topology::FieldBase::SCALAR;
  const int fiberDimE = 1;
  const PylithScalar fieldValuesE[] = {
    sqrt(pow(1.1, 2) + pow(1.2, 2)),
    sqrt(pow(2.1, 2) + pow(2.2, 2)),
    sqrt(pow(3.1, 2) + pow(3.2, 2)),
    sqrt(pow(4.1, 2) + pow(4.2, 2))
  };

  MeshIOAscii iohandler;
  topology::Mesh mesh;
  iohandler.filename(filename);
  iohandler.read(&mesh);

  // Set vertex field
  MeshField field(mesh);
  field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();
  field.vectorFieldType(fieldType);
  field.label(label.c_str());

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  
  const ALE::Obj<RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  CPPUNIT_ASSERT_EQUAL(nvertices, int(vertices->size()));
  int ipt = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt) {
    const PylithScalar* values = &fieldValues[ipt*fiberDim];
    section->updatePoint(*v_iter, values);
  } // for

  VertexFilterVecNorm<MeshField> filter;
  const MeshField& fieldF = filter.filter(field);
  const ALE::Obj<RealSection>& sectionF = fieldF.section();
  CPPUNIT_ASSERT(!sectionF.isNull());

  CPPUNIT_ASSERT_EQUAL(fieldTypeE, fieldF.vectorFieldType());
  CPPUNIT_ASSERT_EQUAL(label, std::string(fieldF.label()));

  ipt = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt) {
    CPPUNIT_ASSERT_EQUAL(fiberDimE, sectionF->getFiberDimension(*v_iter));
    const PylithScalar* values = sectionF->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != values);
    const PylithScalar tolerance = 1.0e-06;
    for (int i=0; i < fiberDimE; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				   values[i]/fieldValuesE[ipt*fiberDimE+i],
				   tolerance);
  } // for
} // testFilter


// End of file 
