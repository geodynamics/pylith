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

#include "TestFieldSubMesh.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SubMesh.hh" // USES SubMesh

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldSubMesh );

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
namespace pylith {
  namespace topology {
    namespace _TestFieldSubMesh {
      const int cellDim = 2;
      const int nvertices = 4;
      const int ncells = 2;
      const int ncorners = 3;
      const int cells[] = {
	0, 1, 3,
	0, 3, 2,
      };
      const PylithScalar coordinates[] = {
	0.0, 0.0,
	1.0, 0.0,
	0.0, 1.0,
	1.0, 1.0,
      };
      const char* label = "bc";
      const int groupSize = 3;
      const int groupVertices[] = {
	1, 2, 3
      };
      const int submeshVertices[] = {
	3, 4, 5
      };
    } // _TestFieldSubMesh
  } // topology
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldSubMesh::testConstructor(void)
{ // testConstructor
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);
} // testConstructor

// ----------------------------------------------------------------------
// Test newSection().
void
pylith::topology::TestFieldSubMesh::testSection(void)
{ // testSection
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);

  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(section.isNull());
} // testSection

// ----------------------------------------------------------------------
// Test mesh().
void
pylith::topology::TestFieldSubMesh::testMesh(void)
{ // testMesh
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);

  const SubMesh& mesh2 = field.mesh();
  CPPUNIT_ASSERT_EQUAL(_TestFieldSubMesh::cellDim-1, mesh2.dimension());  
} // testMesh

// ----------------------------------------------------------------------
// Test spaceDim().
void
pylith::topology::TestFieldSubMesh::testSpaceDim(void)
{ // testSpaceDim
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);

  CPPUNIT_ASSERT_EQUAL(_TestFieldSubMesh::cellDim, field.spaceDim());
} // testSpaceDim

// ----------------------------------------------------------------------
// Test newSection().
void
pylith::topology::TestFieldSubMesh::testNewSection(void)
{ // testNewSection
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);

  field.newSection();
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
} // testNewSection

// ----------------------------------------------------------------------
// Test newSection(points).
void
pylith::topology::TestFieldSubMesh::testNewSectionPoints(void)
{ // testNewSectionPoints
  const int fiberDim = 2;

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  const ALE::Obj<SubMesh::SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  Field<SubMesh> field(submesh);
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  field.newSection(vertices, fiberDim);
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
} // testNewSectionPoints

// ----------------------------------------------------------------------
// Test newSection(domain).
void
pylith::topology::TestFieldSubMesh::testNewSectionDomain(void)
{ // testNewSectionDomain
  const int fiberDim = 2;

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  const ALE::Obj<SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);

  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
} // testNewSectionDomain

// ----------------------------------------------------------------------
// Test newSection(field).
void
pylith::topology::TestFieldSubMesh::testNewSectionField(void)
{ // testNewSectionField
  const int fiberDim = 3;
    
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  const ALE::Obj<SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  // Create field with atlas to use to create new field
  Field<SubMesh> fieldSrc(submesh);
  fieldSrc.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  const ALE::Obj<SubMesh::RealSection>& sectionSrc = fieldSrc.section();
  CPPUNIT_ASSERT(!sectionSrc.isNull());

  const int fiberDim2 = 4;
  Field<SubMesh> field(submesh);
  field.newSection(fieldSrc, fiberDim2);
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim2, section->getFiberDimension(*v_iter));
} // testNewSectionChart

// ----------------------------------------------------------------------
// Test cloneSection().
void
pylith::topology::TestFieldSubMesh::testCloneSection(void)
{ // testCloneSection
  const int fiberDim = 3;
  const int nconstraints[] = { 0, 2, 1, 3 };
  const int constraints[] = {
              // 0
    0, 2,     // 1
    2,        // 2
    0, 1, 2,  // 3
  };
    
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  const ALE::Obj<SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  // Create field with atlas to use to create new field
  Field<SubMesh> fieldSrc(submesh);
  { // Setup source field
    fieldSrc.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
    const ALE::Obj<SubMesh::RealSection>& section = fieldSrc.section();
    CPPUNIT_ASSERT(!section.isNull());
    int iV=0;

    CPPUNIT_ASSERT(!vertices.isNull());
    for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter)
      section->addConstraintDimension(*v_iter, nconstraints[iV++]);
    fieldSrc.allocate();

    int index = 0;
    int i = 0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter, index += nconstraints[i++])
      section->setConstraintDof(*v_iter, &constraints[index]);
    fieldSrc.zero();
  } // Setup source field


  Field<SubMesh> field(submesh);
  field.cloneSection(fieldSrc);
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  int iV = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
    CPPUNIT_ASSERT_EQUAL(nconstraints[iV++], 
			 section->getConstraintDimension(*v_iter));
  } // for
} // testCloneSection

// ----------------------------------------------------------------------
// Test clear().
void
pylith::topology::TestFieldSubMesh::testClear(void)
{ // testClear
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);

  field.scale(2.0);
  field.vectorFieldType(Field<SubMesh>::TENSOR);
  field.addDimensionOkay(true);
  
  field.clear();

  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0), field._metadata.scale);
  CPPUNIT_ASSERT_EQUAL(Field<SubMesh>::OTHER, field._metadata.vectorFieldType);
  CPPUNIT_ASSERT_EQUAL(false, field._metadata.dimsOkay);
} // testClear

// ----------------------------------------------------------------------
// Test allocate().
void
pylith::topology::TestFieldSubMesh::testAllocate(void)
{ // testAllocate
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  const ALE::Obj<SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  scalar_array values(fiberDim);
  int i = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  const PylithScalar tolerance = 1.0e-6;
  i = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], values[iDim], tolerance);
    } // for
  } // for
} // testAllocate

// ----------------------------------------------------------------------
// Test zero().
void
pylith::topology::TestFieldSubMesh::testZero(void)
{ // testZero
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  const ALE::Obj<SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  scalar_array values(fiberDim);
  int i = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.zero();

  const PylithScalar tolerance = 1.0e-6;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[iDim], tolerance);
    } // for
  } // for
} // testZero

// ----------------------------------------------------------------------
// Test complete().
void
pylith::topology::TestFieldSubMesh::testComplete(void)
{ // testComplete
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  const ALE::Obj<SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  scalar_array values(fiberDim);
  int i = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.complete();

  // Expect no change for this serial test
  i = 0;
  const PylithScalar tolerance = 1.0e-6;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], values[iDim], tolerance);
    } // for
  } // for
} // testComplete

// ----------------------------------------------------------------------
// Test copy().
void
pylith::topology::TestFieldSubMesh::testCopy(void)
{ // testCopy
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  const ALE::Obj<SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  Field<SubMesh> fieldSrc(submesh);
  { // Setup source field
    fieldSrc.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    const ALE::Obj<SubMesh::RealSection>& section = fieldSrc.section();
    CPPUNIT_ASSERT(!section.isNull());
    
    scalar_array values(fiberDim);
    int i = 0;
    for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter) {
      for (int iDim=0; iDim < fiberDim; ++iDim)
	values[iDim] = valuesNondim[i++];
      section->updatePoint(*v_iter, &values[0]);
    } // for
  } // Setup source field

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  field.copy(fieldSrc);

  int i = 0;
  scalar_array values(fiberDim);
  const PylithScalar tolerance = 1.0e-6;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], values[iDim], tolerance);
    } // for
  } // for
} // testCopy

// ----------------------------------------------------------------------
// Test operator+=().
void
pylith::topology::TestFieldSubMesh::testOperatorAdd(void)
{ // testOperateAdd
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesA[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };
  const PylithScalar valuesB[] = {
    10.1, 20.2, 30.3,
    10.2, 20.3, 30.4,
    10.3, 20.4, 30.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  const ALE::Obj<SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  Field<SubMesh> fieldSrc(submesh);
  { // Setup source field
    fieldSrc.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    const ALE::Obj<SubMesh::RealSection>& section = fieldSrc.section();
    CPPUNIT_ASSERT(!section.isNull());
    
    scalar_array values(fiberDim);
    int i = 0;
    for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter) {
      for (int iDim=0; iDim < fiberDim; ++iDim)
	values[iDim] = valuesA[i++];
      section->updatePoint(*v_iter, &values[0]);
    } // for
  } // Setup source field

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  { // Setup destination field

    scalar_array values(fiberDim);
    int i = 0;
    for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter) {
      for (int iDim=0; iDim < fiberDim; ++iDim)
	values[iDim] = valuesB[i++];
      section->updatePoint(*v_iter, &values[0]);
    } // for
  } // Setup destination field

  field += fieldSrc;

  int i = 0;
  scalar_array values(fiberDim);
  const PylithScalar tolerance = 1.0e-6;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const PylithScalar valueE = valuesA[i] + valuesB[i];
      ++i;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, values[iDim], tolerance);
    } // for
  } // for
} // testOperateAdd

// ----------------------------------------------------------------------
// Test dimensionalize().
void
pylith::topology::TestFieldSubMesh::testDimensionalize(void)
{ // testDimensionalize
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  const ALE::Obj<SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  scalar_array values(fiberDim);
  int i = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.scale(scale);
  field.addDimensionOkay(true);
  field.dimensionalize();

  i = 0;
  const PylithScalar tolerance = 1.0e-6;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const PylithScalar valueE = valuesNondim[i++]*scale;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, values[iDim], tolerance);
    } // for
  } // for

} // testDimensionalize

// ----------------------------------------------------------------------
// Test view().
void
pylith::topology::TestFieldSubMesh::testView(void)
{ // testView
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  const ALE::Obj<SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  scalar_array values(fiberDim);
  int i = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.view("Testing view");
} // testView

// ----------------------------------------------------------------------
// Test createScatter().
void
pylith::topology::TestFieldSubMesh::testCreateScatter(void)
{ // testCreateScatter
  const int fiberDim = 3;

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);
  field.newSection(FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  const ALE::Obj<SubMesh::SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SubMesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const int sizeE = vertices->size() * fiberDim;

  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatter(submesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field<SubMesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.scatter);
  CPPUNIT_ASSERT(sinfo.scatterVec);
  CPPUNIT_ASSERT(sinfo.vector);

  int size = 0;
  VecGetSize(sinfo.vector, &size);
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  // Make sure we can do multiple calls to createScatter().
  field.createScatter(submesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  // Create another scatter.
  field.createScatter(submesh, "B");
  CPPUNIT_ASSERT_EQUAL(size_t(2), field._scatters.size());
  const Field<SubMesh>::ScatterInfo& sinfoB = field._getScatter("B");
  CPPUNIT_ASSERT(sinfoB.scatter);
  CPPUNIT_ASSERT(sinfoB.scatterVec);
  CPPUNIT_ASSERT(sinfoB.vector);

  Field<SubMesh> field2(submesh);
  field2.cloneSection(field);
  CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

  const Field<SubMesh>::ScatterInfo& sinfo2 = field2._getScatter("");
  CPPUNIT_ASSERT(sinfo2.scatter);
  CPPUNIT_ASSERT(sinfo2.scatterVec);
  CPPUNIT_ASSERT(sinfo2.vector);

  const Field<SubMesh>::ScatterInfo& sinfo2B = field2._getScatter("B");
  CPPUNIT_ASSERT(sinfo2B.scatter);
  CPPUNIT_ASSERT(sinfo2B.scatterVec);
  CPPUNIT_ASSERT(sinfo2B.vector);
} // testCreateScatter

// ----------------------------------------------------------------------
// Test createScatterWithBC().
void
pylith::topology::TestFieldSubMesh::testCreateScatterWithBC(void)
{ // testCreateScatterWithBC
  const int fiberDim = 3;

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);
  field.newSection(FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  const ALE::Obj<SubMesh::SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SubMesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const int sizeE = vertices->size() * fiberDim;

  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatterWithBC(submesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field<SubMesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.scatter);
  CPPUNIT_ASSERT(sinfo.scatterVec);
  CPPUNIT_ASSERT(sinfo.vector);

  int size = 0;
  VecGetSize(sinfo.vector, &size);
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  // Make sure we can do multiple calls to createScatterWithBC().
  field.createScatterWithBC(submesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  // Create another scatter.
  field.createScatterWithBC(submesh, "B");
  CPPUNIT_ASSERT_EQUAL(size_t(2), field._scatters.size());
  const Field<SubMesh>::ScatterInfo& sinfoB = field._getScatter("B");
  CPPUNIT_ASSERT(sinfoB.scatter);
  CPPUNIT_ASSERT(sinfoB.scatterVec);
  CPPUNIT_ASSERT(sinfoB.vector);

  Field<SubMesh> field2(submesh);
  field2.cloneSection(field);
  CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

  const Field<SubMesh>::ScatterInfo& sinfo2 = field2._getScatter("");
  CPPUNIT_ASSERT(sinfo2.scatter);
  CPPUNIT_ASSERT(sinfo2.scatterVec);
  CPPUNIT_ASSERT(sinfo2.vector);

  const Field<SubMesh>::ScatterInfo& sinfo2B = field2._getScatter("B");
  CPPUNIT_ASSERT(sinfo2B.scatter);
  CPPUNIT_ASSERT(sinfo2B.scatterVec);
  CPPUNIT_ASSERT(sinfo2B.vector);
} // testCreateScatterWithBC

// ----------------------------------------------------------------------
// Test vector().
void
pylith::topology::TestFieldSubMesh::testVector(void)
{ // testVector
  const int fiberDim = 3;

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  Field<SubMesh> field(submesh);
  field.newSection(FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();

  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatter(submesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field<SubMesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.scatter);
  CPPUNIT_ASSERT(sinfo.scatterVec);
  CPPUNIT_ASSERT(sinfo.vector);

  const ALE::Obj<SubMesh::SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SubMesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  
  const PetscVec vec = field.vector();
  CPPUNIT_ASSERT_EQUAL(sinfo.vector, vec);
  int size = 0;
  VecGetSize(vec, &size);
  const int sizeE = vertices->size() * fiberDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);
} // testVector

// ----------------------------------------------------------------------
// Test scatterSectionToVector().
void
pylith::topology::TestFieldSubMesh::testScatterSectionToVector(void)
{ // testScatterSectionToVector
  const char* context = "abc";
  const int fiberDim = 3;
  const PylithScalar valuesE[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);
  const ALE::Obj<SubMesh::SieveMesh>& sieveMesh = submesh.sieveMesh();
  field.newSection(FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  const ALE::Obj<SubMesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  scalar_array values(fiberDim);
  int i = 0;
  for (SubMesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesE[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.createScatter(submesh, context);
  field.scatterSectionToVector(context);
  const PetscVec vec = field.vector(context);
  CPPUNIT_ASSERT(0 != vec);
  int size = 0;
  VecGetSize(vec, &size);
  PylithScalar* valuesVec = 0;
  VecGetArray(vec, &valuesVec);

  const PylithScalar tolerance = 1.0e-06;
  const int sizeE = vertices->size() * fiberDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);
  for (int i=0; i < sizeE; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i], valuesVec[i], tolerance);
  VecRestoreArray(vec, &valuesVec);
} // testScatterSectionToVector

// ----------------------------------------------------------------------
// Test scatterVectorToSection().
void
pylith::topology::TestFieldSubMesh::testScatterVectorToSection(void)
{ // testScatterVectorToSection
  const char* context = "abcd";
  const int fiberDim = 3;
  const PylithScalar valuesE[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);
  field.newSection(FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();
  field.createScatter(submesh, context);

  const ALE::Obj<SubMesh::SieveMesh>& sieveMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SubMesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  const PetscVec vec = field.vector(context);
  CPPUNIT_ASSERT(0 != vec);
  int size = 0;
  VecGetSize(vec, &size);
  const int sizeE = vertices->size() * fiberDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  const PylithScalar tolerance = 1.0e-06;
  PylithScalar* valuesVec = 0;
  VecGetArray(vec, &valuesVec);
  for (int i=0; i < sizeE; ++i)
    valuesVec[i] = valuesE[i];
  VecRestoreArray(vec, &valuesVec);

  field.scatterVectorToSection(context);

  scalar_array values(fiberDim);
  int i = 0;
  const ALE::Obj<SubMesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  for (SubMesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], fiberDim);
    for (int iDim=0; iDim < fiberDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i++], values[iDim], tolerance);
  } // for
} // testScatterVectorToSection

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldSubMesh::_buildMesh(Mesh* mesh,
					       SubMesh* submesh)
{ // _buildMesh
  assert(0 != mesh);
  assert(0 != submesh);

  mesh->createSieveMesh(_TestFieldSubMesh::cellDim);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  ALE::Obj<Mesh::SieveMesh::sieve_type> sieve = 
    new Mesh::SieveMesh::sieve_type(sieveMesh->comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  ALE::Obj<SieveFlexMesh::sieve_type> s = 
    new SieveFlexMesh::sieve_type(sieve->comm(), sieve->debug());
  CPPUNIT_ASSERT(!s.isNull());
  
  const int cellDim = _TestFieldSubMesh::cellDim;
  const int ncells = _TestFieldSubMesh::ncells;
  const int* cells = _TestFieldSubMesh::cells;
  const int nvertices = _TestFieldSubMesh::nvertices;
  const int ncorners = _TestFieldSubMesh::ncorners;
  const int spaceDim = _TestFieldSubMesh::cellDim;
  const PylithScalar* coordinates = _TestFieldSubMesh::coordinates;
  const bool interpolate = false;
  ALE::SieveBuilder<SieveFlexMesh>::buildTopology(s, cellDim, ncells, (int*) cells,
					      nvertices, interpolate, 
					      ncorners);
  std::map<Mesh::SieveMesh::point_type,Mesh::SieveMesh::point_type> renumbering;
  ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
  sieveMesh->setSieve(sieve);
  sieveMesh->stratify();
  ALE::SieveBuilder<SieveMesh>::buildCoordinates(sieveMesh, spaceDim, 
						 coordinates);

  typedef Mesh::SieveMesh::int_section_type::chart_type chart_type;
  const ALE::Obj<SieveMesh::int_section_type>& groupField = 
    sieveMesh->getIntSection(_TestFieldSubMesh::label);
  assert(!groupField.isNull());

  const int numPoints = _TestFieldSubMesh::groupSize;
  const int numVertices = sieveMesh->depthStratum(0)->size();
  const int numCells = sieveMesh->heightStratum(0)->size();
  groupField->setChart(chart_type(numCells, numCells+numVertices));
  for(int i=0; i < numPoints; ++i)
    groupField->setFiberDimension(numCells+_TestFieldSubMesh::groupVertices[i],
				  1);
  sieveMesh->allocate(groupField);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);

  submesh->createSubMesh(*mesh, _TestFieldSubMesh::label);
} // _buildMesh


// End of file 
