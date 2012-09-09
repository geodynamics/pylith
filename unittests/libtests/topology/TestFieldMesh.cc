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

#include "TestFieldMesh.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldMesh );

// ----------------------------------------------------------------------
namespace pylith {
  namespace topology {
    namespace _TestFieldMesh {
      const int cellDim = 2;
      const int nvertices = 4;
      const int ncells = 1;
      const int ncorners = 4;
      const int cells[] = { 0, 1, 2, 3 };
      const PylithScalar coordinates[] = {
	0.0, 0.0,
	1.0, 0.0,
	0.0, 1.0,
	1.0, 1.0,
      };
    } // _TestFieldMesh
  } // topology
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldMesh::testConstructor(void)
{ // testConstructor
  Mesh mesh;
  Field<Mesh> field(mesh);
} // testConstructor

// ----------------------------------------------------------------------
// Test section().
void
pylith::topology::TestFieldMesh::testSection(void)
{ // testSection
  Mesh mesh;
  Field<Mesh> field(mesh);

  mesh.createSieveMesh();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(section.isNull());
} // testSection

// ----------------------------------------------------------------------
// Test mesh().
void
pylith::topology::TestFieldMesh::testMesh(void)
{ // testMesh
  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);

  const Mesh& mesh2 = field.mesh();
  CPPUNIT_ASSERT_EQUAL(_TestFieldMesh::cellDim, mesh2.dimension());  
} // testMesh

// ----------------------------------------------------------------------
// Test label().
void
pylith::topology::TestFieldMesh::testLabel(void)
{ // testLabel
  const std::string label = "velocity";

  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);

  field.label(label.c_str());
  CPPUNIT_ASSERT_EQUAL(label, std::string(field.label()));
} // testLabel

// ----------------------------------------------------------------------
// Test vectorFieldType().
void
pylith::topology::TestFieldMesh::testVectorFieldType(void)
{ // testVectorFieldType
} // testVectorFieldType

// ----------------------------------------------------------------------
// Test scale().
void
pylith::topology::TestFieldMesh::testScale(void)
{ // testScale
} // testScale

// ----------------------------------------------------------------------
// Test addDimensionsOkay().
void
pylith::topology::TestFieldMesh::testAddDimensionsOkay(void)
{ // testAddDimensionsOkay
} // testAddDimensionsOkay

// ----------------------------------------------------------------------
// Test spaceDim().
void
pylith::topology::TestFieldMesh::testSpaceDim(void)
{ // testSpaceDim
  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);

  CPPUNIT_ASSERT_EQUAL(_TestFieldMesh::cellDim, field.spaceDim());
} // testSpaceDim

// ----------------------------------------------------------------------
// Test newSection().
void
pylith::topology::TestFieldMesh::testNewSection(void)
{ // testNewSection
  Mesh mesh;
  _buildMesh(&mesh);

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());

  field.newSection();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  CPPUNIT_ASSERT_EQUAL(label, std::string(section->getName()));
} // testNewSection

// ----------------------------------------------------------------------
// Test newSection(points).
void
pylith::topology::TestFieldMesh::testNewSectionPoints(void)
{ // testNewSectionPoints
  const int fiberDim = 2;

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());

  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  field.newSection(vertices, fiberDim);
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  CPPUNIT_ASSERT(!vertices.isNull());
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));

  CPPUNIT_ASSERT_EQUAL(label, std::string(section->getName()));
} // testNewSectionPoints

// ----------------------------------------------------------------------
// Test newSection(int_array).
void
pylith::topology::TestFieldMesh::testNewSectionPointsArray(void)
{ // testNewSectionPointsArray
  const int fiberDim = 2;

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());

  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const int npts = vertices->size() / 2;
  int_array pointsIn(npts);
  int_array pointsOut(vertices->size() - npts);
  int count = 0;
  size_t iIn = 0;
  size_t iOut = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    if (count % 2  == 0)
      pointsIn[iIn++] = *v_iter;
    else
      pointsOut[iOut++] = *v_iter;
    ++count;
  } // for
  CPPUNIT_ASSERT_EQUAL(iIn, pointsIn.size());
  CPPUNIT_ASSERT_EQUAL(iOut, pointsOut.size());

  field.newSection(pointsIn, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  // Points in array should have a fiber dimension of fiberDim.
  for (int i=0; i < pointsIn.size(); ++i)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(pointsIn[i]));
  
  // Points not int array should have a fiber dimension of zero.
  for (int i=0; i < pointsOut.size(); ++i)
    CPPUNIT_ASSERT_EQUAL(0, section->getFiberDimension(pointsOut[i]));

  CPPUNIT_ASSERT_EQUAL(label, std::string(section->getName()));
} // testNewSectionPointsArray

// ----------------------------------------------------------------------
// Test newSection(domain).
void
pylith::topology::TestFieldMesh::testNewSectionDomain(void)
{ // testNewSectionDomain
  const int fiberDim = 2;

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);

  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));

  CPPUNIT_ASSERT_EQUAL(label, std::string(section->getName()));
} // testNewSectionDomain

// ----------------------------------------------------------------------
// Test newSection(field).
void
pylith::topology::TestFieldMesh::testNewSectionField(void)
{ // testNewSectionField
  const int fiberDim = 3;
    
  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  const ALE::Obj<Mesh::RealSection>& sectionSrc = fieldSrc.section();
  CPPUNIT_ASSERT(!sectionSrc.isNull());
  const Mesh::RealSection::chart_type& chart = sectionSrc->getChart();

  const int fiberDim2 = 5;
  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(fieldSrc, fiberDim2);
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim2, section->getFiberDimension(*v_iter));

  CPPUNIT_ASSERT_EQUAL(label, std::string(section->getName()));
} // testNewSectionField

// ----------------------------------------------------------------------
// Test cloneSection().
void
pylith::topology::TestFieldMesh::testCloneSection(void)
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
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
    CPPUNIT_ASSERT(!section.isNull());
    int iV=0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
    fieldSrc.createScatter(mesh);
    fieldSrc.createScatter(mesh, "A");
  } // Setup source field

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.cloneSection(fieldSrc);
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  int iV = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
    CPPUNIT_ASSERT_EQUAL(nconstraints[iV++], 
			 section->getConstraintDimension(*v_iter));
  } // for

  // Verify vector scatters were also copied.
  CPPUNIT_ASSERT_EQUAL(fieldSrc._scatters[""].scatter, 
		       field._scatters[""].scatter);
  CPPUNIT_ASSERT_EQUAL(fieldSrc._scatters["A"].scatter, 
		       field._scatters["A"].scatter);

  CPPUNIT_ASSERT_EQUAL(label, std::string(section->getName()));
} // testCloneSection

// ----------------------------------------------------------------------
// Test clear().
void
pylith::topology::TestFieldMesh::testClear(void)
{ // testClear
  Mesh mesh(_TestFieldMesh::cellDim);
  Field<Mesh> field(mesh);

  field.scale(2.0);
  field.vectorFieldType(Field<Mesh>::TENSOR);
  field.addDimensionOkay(true);
  
  field.clear();

  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0), field._metadata["default"].scale);
  CPPUNIT_ASSERT_EQUAL(Field<Mesh>::OTHER, field._metadata["default"].vectorFieldType);
  CPPUNIT_ASSERT_EQUAL(false, field._metadata["default"].dimsOkay);
} // testClear

// ----------------------------------------------------------------------
// Test allocate().
void
pylith::topology::TestFieldMesh::testAllocate(void)
{ // testAllocate
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  scalar_array values(fiberDim);
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  const PylithScalar tolerance = 1.0e-6;
  i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::topology::TestFieldMesh::testZero(void)
{ // testZero
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  scalar_array values(fiberDim);
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.zero();

  const PylithScalar tolerance = 1.0e-6;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[iDim], tolerance);
    } // for
  } // for
} // testZero

// ----------------------------------------------------------------------
// Test zero().
void
pylith::topology::TestFieldMesh::testZeroAll(void)
{ // testZeroAll
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };
  const int nconstraints[] = { 0, 2, 1, 3 };
  const int constraints[] = {
              // 0
    0, 2,     // 1
    2,        // 2
    0, 1, 2,  // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  // Create field and set constraint sizes
  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  int iV=0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    section->addConstraintDimension(*v_iter, nconstraints[iV++]);
  field.allocate();
  int index = 0;
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter, index += nconstraints[i++])
    section->setConstraintDof(*v_iter, &constraints[index]);
  field.zero();

  scalar_array values(fiberDim);
  i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePointAll(*v_iter, &values[0]);
  } // for
  
  field.zeroAll();
  
  const PylithScalar tolerance = 1.0e-6;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[iDim], tolerance);
    } // for
  } // for
} // testZeroAll

// ----------------------------------------------------------------------
// Test complete().
void
pylith::topology::TestFieldMesh::testComplete(void)
{ // testComplete
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  scalar_array values(fiberDim);
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::topology::TestFieldMesh::testCopy(void)
{ // testCopy
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
    CPPUNIT_ASSERT(!section.isNull());
    
    scalar_array values(fiberDim);
    int i = 0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter) {
      for (int iDim=0; iDim < fiberDim; ++iDim)
	values[iDim] = valuesNondim[i++];
      section->updatePoint(*v_iter, &values[0]);
    } // for
  } // Setup source field

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  field.copy(fieldSrc);

  int i = 0;
  scalar_array values(fiberDim);
  const PylithScalar tolerance = 1.0e-6;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::topology::TestFieldMesh::testOperatorAdd(void)
{ // testOperateAdd
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesA[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };
  const PylithScalar valuesB[] = {
    10.1, 20.2, 30.3,
    10.2, 20.3, 30.4,
    10.3, 20.4, 30.5,
    10.4, 20.5, 30.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
    CPPUNIT_ASSERT(!section.isNull());
    
    scalar_array values(fiberDim);
    int i = 0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter) {
      for (int iDim=0; iDim < fiberDim; ++iDim)
	values[iDim] = valuesA[i++];
      section->updatePoint(*v_iter, &values[0]);
    } // for
  } // Setup source field

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  { // Setup destination field

    scalar_array values(fiberDim);
    int i = 0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::topology::TestFieldMesh::testDimensionalize(void)
{ // testDimensionalize
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  scalar_array values(fiberDim);
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::topology::TestFieldMesh::testView(void)
{ // testView
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  scalar_array values(fiberDim);
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::topology::TestFieldMesh::testCreateScatter(void)
{ // testCreateScatter
  const int fiberDim = 3;

  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const int sizeE = vertices->size() * fiberDim;

  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatter(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field<Mesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.scatter);
  CPPUNIT_ASSERT(sinfo.scatterVec);
  CPPUNIT_ASSERT(sinfo.vector);

  int size = 0;
  VecGetSize(sinfo.vector, &size);
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  // Check vector name
  const char* vecname = 0;
  PetscObjectGetName((PetscObject)sinfo.vector, &vecname);
  CPPUNIT_ASSERT_EQUAL(label, std::string(vecname));

  // Make sure we can do multiple calls to createScatter().
  field.createScatter(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  // Create another scatter.
  field.createScatter(mesh, "B");
  CPPUNIT_ASSERT_EQUAL(size_t(2), field._scatters.size());
  const Field<Mesh>::ScatterInfo& sinfoB = field._getScatter("B");
  CPPUNIT_ASSERT(sinfoB.scatter);
  CPPUNIT_ASSERT(sinfoB.scatterVec);
  CPPUNIT_ASSERT(sinfoB.vector);

  Field<Mesh> field2(mesh);
  field2.cloneSection(field);
  CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

  const Field<Mesh>::ScatterInfo& sinfo2 = field2._getScatter("");
  CPPUNIT_ASSERT(sinfo2.scatter);
  CPPUNIT_ASSERT(sinfo2.scatterVec);
  CPPUNIT_ASSERT(sinfo2.vector);

  const Field<Mesh>::ScatterInfo& sinfo2B = field2._getScatter("B");
  CPPUNIT_ASSERT(sinfo2B.scatter);
  CPPUNIT_ASSERT(sinfo2B.scatterVec);
  CPPUNIT_ASSERT(sinfo2B.vector);

} // testCreateScatter

// ----------------------------------------------------------------------
// Test createScatterWithBC().
void
pylith::topology::TestFieldMesh::testCreateScatterWithBC(void)
{ // testCreateScatterWithBC
  const int fiberDim = 3;

  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const int sizeE = vertices->size() * fiberDim;

  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatterWithBC(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field<Mesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.scatter);
  CPPUNIT_ASSERT(sinfo.scatterVec);
  CPPUNIT_ASSERT(sinfo.vector);

  int size = 0;
  VecGetSize(sinfo.vector, &size);
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  // Check vector name
  const char* vecname = 0;
  PetscObjectGetName((PetscObject)sinfo.vector, &vecname);
  CPPUNIT_ASSERT_EQUAL(label, std::string(vecname));

  // Make sure we can do multiple calls to createScatterWithBC().
  field.createScatterWithBC(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  // Create another scatter.
  field.createScatterWithBC(mesh, "B");
  CPPUNIT_ASSERT_EQUAL(size_t(2), field._scatters.size());
  const Field<Mesh>::ScatterInfo& sinfoB = field._getScatter("B");
  CPPUNIT_ASSERT(sinfoB.scatter);
  CPPUNIT_ASSERT(sinfoB.scatterVec);
  CPPUNIT_ASSERT(sinfoB.vector);

  Field<Mesh> field2(mesh);
  field2.cloneSection(field);
  CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

  const Field<Mesh>::ScatterInfo& sinfo2 = field2._getScatter("");
  CPPUNIT_ASSERT(sinfo2.scatter);
  CPPUNIT_ASSERT(sinfo2.scatterVec);
  CPPUNIT_ASSERT(sinfo2.vector);

  const Field<Mesh>::ScatterInfo& sinfo2B = field2._getScatter("B");
  CPPUNIT_ASSERT(sinfo2B.scatter);
  CPPUNIT_ASSERT(sinfo2B.scatterVec);
  CPPUNIT_ASSERT(sinfo2B.vector);

} // testCreateScatterWithBC

// ----------------------------------------------------------------------
// Test vector().
void
pylith::topology::TestFieldMesh::testVector(void)
{ // testVector
  const int fiberDim = 3;

  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatter(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field<Mesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.scatter);
  CPPUNIT_ASSERT(sinfo.scatterVec);
  CPPUNIT_ASSERT(sinfo.vector);

  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
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
pylith::topology::TestFieldMesh::testScatterSectionToVector(void)
{ // testScatterSectionToVector
  const char* context = "abc";
  const int fiberDim = 3;
  const PylithScalar valuesE[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  scalar_array values(fiberDim);
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesE[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.createScatter(mesh, context);
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
pylith::topology::TestFieldMesh::testScatterVectorToSection(void)
{ // testScatterVectorToSection
  const char* context = "abcd";
  const int fiberDim = 3;
  const PylithScalar valuesE[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  field.createScatter(mesh, context);
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
    valuesVec[i] = valuesE[i];
  VecRestoreArray(vec, &valuesVec);

  field.createScatter(mesh, context);
  field.scatterVectorToSection(context);

  scalar_array values(fiberDim);
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], fiberDim);
    for (int iDim=0; iDim < fiberDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i++], values[iDim], tolerance);
  } // for
} // testScatterVectorToSection

// ----------------------------------------------------------------------
// Test splitDefault().
void
pylith::topology::TestFieldMesh::testSplitDefault(void)
{ // testSplitDefault
  const int spaceDim = _TestFieldMesh::cellDim;
  const int numFibrations = spaceDim;
  const int nconstraints[4] = { 1, 2, 0, 1 };
  const int constraints[4] = {
    1,     // 0
    0, 1,  // 1
           // 2
    0,     // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, spaceDim);
    fieldSrc.splitDefault();
    const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
    CPPUNIT_ASSERT(!section.isNull());
    int iV=0;
    int iC=0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
        v_iter != vertices->end();
        ++v_iter, ++iV) {
      const int nconstraintsVertex = nconstraints[iV];
      section->addConstraintDimension(*v_iter, nconstraintsVertex);
      for (int iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const int fibration = constraints[iC++];
        section->addConstraintDimension(*v_iter, 1, fibration);
      } // for
    } // for
    fieldSrc.allocate();

    iC = 0;
    iV = 0;
    int zero = 0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
        v_iter != vertices->end();
        ++v_iter, ++iV) {
      const int nconstraintsVertex = nconstraints[iV];
      if (nconstraintsVertex > 0)
        section->setConstraintDof(*v_iter, &constraints[iC]);
      for (int iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const int fibration = constraints[iC++];
        section->setConstraintDof(*v_iter, &zero, fibration);
      } // for
    } // for
  } // Setup source field

  const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
  CPPUNIT_ASSERT(!section.isNull());
  CPPUNIT_ASSERT_EQUAL(numFibrations, section->getNumSpaces());

  for (int fibration=0; fibration < spaceDim; ++fibration) {
    const ALE::Obj<Mesh::RealSection>& sectionSplit = section->getFibration(fibration);
    CPPUNIT_ASSERT(!sectionSplit.isNull());
    CPPUNIT_ASSERT(!vertices.isNull());
    int iV = 0;
    int iC = 0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
        v_iter != vertices->end();
        ++v_iter, ++iV) {
      CPPUNIT_ASSERT_EQUAL(1, section->getFiberDimension(*v_iter, fibration));
      bool isConstrained = false;
      const int nconstraintsVertex = nconstraints[iV];
      for (int iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint)
        if (constraints[iC++] == fibration)
          isConstrained = true;
      const int constraintDimE = (!isConstrained) ? 0 : 1;
      CPPUNIT_ASSERT_EQUAL(constraintDimE,
                           section->getConstraintDimension(*v_iter, fibration));
    } // for
  } // for
} // testSplitDefault

// ----------------------------------------------------------------------
// Test cloneSection() with split field.
void
pylith::topology::TestFieldMesh::testCloneSectionSplit(void)
{ // testCloneSectionSplit
  const int spaceDim = _TestFieldMesh::cellDim;
  const int numFibrations = spaceDim;
  const int nconstraints[4] = { 1, 2, 0, 1 };
  const int constraints[4] = {
    1,     // 0
    0, 1,  // 1
           // 2
    0,     // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, spaceDim);
    fieldSrc.splitDefault();
    const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
    CPPUNIT_ASSERT(!section.isNull());
    int iV=0;
    int iC=0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
        v_iter != vertices->end();
        ++v_iter, ++iV) {
      const int nconstraintsVertex = nconstraints[iV];
      section->addConstraintDimension(*v_iter, nconstraintsVertex);
      for (int iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const int fibration = constraints[iC++];
        section->addConstraintDimension(*v_iter, 1, fibration);
      } // for
    } // for
    fieldSrc.allocate();

    iC = 0;
    iV = 0;
    int zero = 0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
        v_iter != vertices->end();
        ++v_iter, ++iV) {
      const int nconstraintsVertex = nconstraints[iV];
      if (nconstraintsVertex > 0)
        section->setConstraintDof(*v_iter, &constraints[iC]);
      for (int iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const int fibration = constraints[iC++];
        section->setConstraintDof(*v_iter, &zero, fibration);
      } // for
    } // for
  } // Setup source field

  Field<Mesh> field(mesh);
  field.cloneSection(fieldSrc);

  const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
  CPPUNIT_ASSERT(!section.isNull());
  CPPUNIT_ASSERT_EQUAL(numFibrations, section->getNumSpaces());

  for (int fibration=0; fibration < spaceDim; ++fibration) {
    const ALE::Obj<Mesh::RealSection>& sectionSplit = section->getFibration(fibration);
    CPPUNIT_ASSERT(!sectionSplit.isNull());
    CPPUNIT_ASSERT(!vertices.isNull());
    int iV = 0;
    int iC = 0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
        v_iter != vertices->end();
        ++v_iter, ++iV) {
      CPPUNIT_ASSERT_EQUAL(1, section->getFiberDimension(*v_iter, fibration));
      bool isConstrained = false;
      const int nconstraintsVertex = nconstraints[iV];
      for (int iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint)
        if (constraints[iC++] == fibration)
          isConstrained = true;
      const int constraintDimE = (!isConstrained) ? 0 : 1;
      CPPUNIT_ASSERT_EQUAL(constraintDimE,
                           section->getConstraintDimension(*v_iter, fibration));
    } // for
  } // for
} // testCloneSectionSplit

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_buildMesh(Mesh* mesh)
{ // _buildMesh
  assert(0 != mesh);

  mesh->createSieveMesh(_TestFieldMesh::cellDim);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh->sieveMesh();

  ALE::Obj<Mesh::SieveMesh::sieve_type> sieve = 
    new Mesh::SieveMesh::sieve_type(sieveMesh->comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  ALE::Obj<SieveFlexMesh::sieve_type> s = 
    new SieveFlexMesh::sieve_type(sieve->comm(), sieve->debug());
  
  const int cellDim = _TestFieldMesh::cellDim;
  const int ncells = _TestFieldMesh::ncells;
  const int* cells = _TestFieldMesh::cells;
  const int nvertices = _TestFieldMesh::nvertices;
  const int ncorners = _TestFieldMesh::ncorners;
  const int spaceDim = _TestFieldMesh::cellDim;
  const PylithScalar* coordinates = _TestFieldMesh::coordinates;
  const bool interpolate = false;
  ALE::SieveBuilder<SieveFlexMesh>::buildTopology(s, cellDim, ncells, (int*) cells,
					      nvertices, interpolate, 
					      ncorners);
  std::map<Mesh::SieveMesh::point_type,Mesh::SieveMesh::point_type> renumbering;
  ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
  sieveMesh->setSieve(sieve);
  sieveMesh->stratify();
  ALE::SieveBuilder<Mesh::SieveMesh>::buildCoordinates(sieveMesh, spaceDim, 
						       coordinates);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);
} // _buildMesh


// End of file 
