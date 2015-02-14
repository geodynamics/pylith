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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFieldMesh.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::createDMMesh()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

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
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  Field field(mesh);

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test section().
void
pylith::topology::TestFieldMesh::testSection(void)
{ // testSection
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  Field field(mesh);

  PetscSection section = field.localSection();
  CPPUNIT_ASSERT(!section);

  PYLITH_METHOD_END;
} // testSection

// ----------------------------------------------------------------------
// Test mesh().
void
pylith::topology::TestFieldMesh::testMesh(void)
{ // testMesh
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _buildMesh(&mesh);
  Field field(mesh);

  const Mesh& mesh2 = field.mesh();
  CPPUNIT_ASSERT_EQUAL(_TestFieldMesh::cellDim, mesh2.dimension());  

  PYLITH_METHOD_END;
} // testMesh

// ----------------------------------------------------------------------
// Test label().
void
pylith::topology::TestFieldMesh::testLabel(void)
{ // testLabel
  PYLITH_METHOD_BEGIN;

  const std::string label = "velocity";

  Mesh mesh;
  _buildMesh(&mesh);
  Field field(mesh);

  field.label(label.c_str());
  CPPUNIT_ASSERT_EQUAL(label, std::string(field.label()));

  PYLITH_METHOD_END;
} // testLabel

// ----------------------------------------------------------------------
// Test vectorFieldType().
void
pylith::topology::TestFieldMesh::testVectorFieldType(void)
{ // testVectorFieldType
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _buildMesh(&mesh);
  Field field(mesh);

  field.vectorFieldType(FieldBase::SCALAR);
  CPPUNIT_ASSERT_EQUAL(FieldBase::SCALAR, field._metadata.vectorFieldType);

  PYLITH_METHOD_END;
} // testVectorFieldType

// ----------------------------------------------------------------------
// Test scale().
void
pylith::topology::TestFieldMesh::testScale(void)
{ // testScale
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _buildMesh(&mesh);
  Field field(mesh);

  const PylithScalar scale = 2.0;
  field.scale(scale);
  const PylithScalar tolerance = 1.0e-6;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(scale, field._metadata.scale, tolerance);

  PYLITH_METHOD_END;
} // testScale

// ----------------------------------------------------------------------
// Test dimensionalizeOkay().
void
pylith::topology::TestFieldMesh::testAddDimensionOkay(void)
{ // testAddDimensionOkay
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _buildMesh(&mesh);
  Field field(mesh);

  CPPUNIT_ASSERT_EQUAL(false, field._metadata.dimsOkay);
  field.dimensionalizeOkay(true);
  CPPUNIT_ASSERT_EQUAL(true, field._metadata.dimsOkay);

  PYLITH_METHOD_END;
} // testAddDimensionOkay

// ----------------------------------------------------------------------
// Test spaceDim().
void
pylith::topology::TestFieldMesh::testSpaceDim(void)
{ // testSpaceDim
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _buildMesh(&mesh);
  Field field(mesh);

  CPPUNIT_ASSERT_EQUAL(_TestFieldMesh::cellDim, field.spaceDim());

  PYLITH_METHOD_END;
} // testSpaceDim

// ----------------------------------------------------------------------
// Test chartSize().
void 
pylith::topology::TestFieldMesh::testChartSize(void)
{ // testChartSize
  PYLITH_METHOD_BEGIN;

  const int fiberDim = 2;
  const std::string& label = "field A";

  Mesh mesh;
  _buildMesh(&mesh);

  Field field(mesh);
  field.label(label.c_str());

  CPPUNIT_ASSERT_EQUAL(0, field.chartSize());

  field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();

  CPPUNIT_ASSERT_EQUAL(_TestFieldMesh::nvertices, field.chartSize());

  PYLITH_METHOD_END;
} // testChartSize

// ----------------------------------------------------------------------
// Test sectionSize().
void 
pylith::topology::TestFieldMesh::testSectionSize(void)
{ // testSectionSize
  PYLITH_METHOD_BEGIN;

  const int fiberDim = 2;
  const std::string& label = "field A";

  Mesh mesh;
  _buildMesh(&mesh);

  Field field(mesh);
  field.label(label.c_str());

  CPPUNIT_ASSERT_EQUAL(0, field.sectionSize());

  field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();

  CPPUNIT_ASSERT_EQUAL(_TestFieldMesh::nvertices*fiberDim, field.sectionSize());

  PYLITH_METHOD_END;
} // testSectionSize

// ----------------------------------------------------------------------
// Test hasSection().
void 
pylith::topology::TestFieldMesh::testHasSection(void)
{ // testHasSection
  PYLITH_METHOD_BEGIN;

  const int fiberDim = 2;
  const std::string& label = "field A";

  Mesh mesh;
  _buildMesh(&mesh);
  Field field(mesh);
  field.label(label.c_str());

  CPPUNIT_ASSERT(!field.hasSection());
  field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
  
  CPPUNIT_ASSERT(field.hasSection());

  PYLITH_METHOD_END;
} // testHasSection

// ----------------------------------------------------------------------
// Test newSection(points).
void
pylith::topology::TestFieldMesh::testNewSectionPoints(void)
{ // testNewSectionPoints
  PYLITH_METHOD_BEGIN;

  const int fiberDim = 2;
  const std::string& label = "field A";

  Mesh mesh;
  _buildMesh(&mesh);

  Field field(mesh);
  field.label(label.c_str());
  field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
  } // for

  const char *name = NULL;
  PetscVec vec = field.localVector();CPPUNIT_ASSERT(vec);
  PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));

  PYLITH_METHOD_END;
} // testNewSectionPoints

// ----------------------------------------------------------------------
// Test newSection(int_array).
void
pylith::topology::TestFieldMesh::testNewSectionPointsArray(void)
{ // testNewSectionPointsArray
  PYLITH_METHOD_BEGIN;

  const int fiberDim = 2;

  Mesh mesh;
  _buildMesh(&mesh);

  Field field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  const int npts = (vEnd-vStart) / 2;
  int_array pointsIn(npts);
  int_array pointsOut(vEnd-vStart - npts);
  int count = 0;
  size_t iIn = 0;
  size_t iOut = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    if (count % 2  == 0)
      pointsIn[iIn++] = v;
    else
      pointsOut[iOut++] = v;
    ++count;
  } // for
  CPPUNIT_ASSERT_EQUAL(iIn, pointsIn.size());
  CPPUNIT_ASSERT_EQUAL(iOut, pointsOut.size());

  field.newSection(pointsIn, fiberDim);
  field.allocate();

  // Points in array should have a fiber dimension of fiberDim.
  VecVisitorMesh fieldVisitor(field);
  for(int i = 0; i < pointsIn.size(); ++i) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(pointsIn[i]));
  } // for

  // Points not int array should have a fiber dimension of zero.
  PetscInt pStart, pEnd;
  PetscSection section = field.localSection();CPPUNIT_ASSERT(section);
  PetscErrorCode err = PetscSectionGetChart(section, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
  for(int i = 0; i < pointsOut.size(); ++i) {
    if (pointsOut[i] >= pStart && pointsOut[i] < pEnd) {
      CPPUNIT_ASSERT_EQUAL(0, fieldVisitor.sectionDof(pointsOut[i]));
    } // if
  } // for

  const char *name = NULL;
  PetscVec vec = field.localVector();CPPUNIT_ASSERT(vec);
  err = PetscObjectGetName((PetscObject) vec, &name);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));

  PYLITH_METHOD_END;
} // testNewSectionPointsArray

// ----------------------------------------------------------------------
// Test newSection(domain).
void
pylith::topology::TestFieldMesh::testNewSectionDomain(void)
{ // testNewSectionDomain
  PYLITH_METHOD_BEGIN;

  const int fiberDim = 2;

  Mesh mesh;
  _buildMesh(&mesh);

  Field field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(Field::VERTICES_FIELD, fiberDim);
  field.allocate();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
  } // for

  const char *name = NULL;
  PetscVec vec = field.localVector();CPPUNIT_ASSERT(vec);
  PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));

  PYLITH_METHOD_END;
} // testNewSectionDomain

// ----------------------------------------------------------------------
// Test newSection(field).
void
pylith::topology::TestFieldMesh::testNewSectionField(void)
{ // testNewSectionField
  PYLITH_METHOD_BEGIN;

  const int fiberDim = 3;
    
  Mesh mesh;
  _buildMesh(&mesh);

  // Create field with atlas to use to create new field
  Field fieldSrc(mesh);
  fieldSrc.newSection(Field::VERTICES_FIELD, fiberDim);
  fieldSrc.allocate();

  const int fiberDim2 = 5;
  Field field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(fieldSrc, fiberDim2);
  field.allocate();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(fiberDim2, fieldVisitor.sectionDof(v));
  } // for

  const char *name = NULL;
  PetscVec vec = field.localVector();CPPUNIT_ASSERT(vec);
  PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));

  PYLITH_METHOD_END;
} // testNewSectionField

// ----------------------------------------------------------------------
// Test cloneSection().
void
pylith::topology::TestFieldMesh::testCloneSection(void)
{ // testCloneSection
  PYLITH_METHOD_BEGIN;

  const PetscInt fiberDim = 3;
  const PetscInt nconstraints[] = { 0, 2, 1, 3 };
  const PetscInt constraints[] = {
              // 0
    0, 2,     // 1
    2,        // 2
    0, 1, 2,  // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  PetscErrorCode err = 0;

  // Create field with atlas to use to create new field
  Field fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field::VERTICES_FIELD, fiberDim);
    PetscSection section = fieldSrc.localSection();CPPUNIT_ASSERT(section);
    for(PetscInt v = vStart, iV = 0; v < vEnd; ++v) {
      err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]);PYLITH_CHECK_ERROR(err);
    } // for
    fieldSrc.allocate();

    int index = 0;
    for(PetscInt v = vStart, iV = 0; v < vEnd; ++v, index += nconstraints[iV++]) {
      err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[index]);PYLITH_CHECK_ERROR(err);
    }
    fieldSrc.zero();
    fieldSrc.createScatter(mesh);
    fieldSrc.createScatter(mesh, "A");
  } // Setup source field

  Field field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.cloneSection(fieldSrc);
  PetscSection section = field.localSection();
  PetscVec vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
  for(PetscInt v = vStart, iV = 0; v < vEnd; ++v) {
    PetscInt dof, cdof;
    err = PetscSectionGetDof(section, v, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetConstraintDof(section, v, &cdof);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDim, dof);
    CPPUNIT_ASSERT_EQUAL(nconstraints[iV++], cdof);
  } // for

  // Verify vector scatters were also copied.
  CPPUNIT_ASSERT_EQUAL(fieldSrc._scatters[""].dm,  field._scatters[""].dm);
  CPPUNIT_ASSERT_EQUAL(fieldSrc._scatters["A"].dm, field._scatters["A"].dm);
  const char *name = NULL;
  err = PetscObjectGetName((PetscObject) vec, &name);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));

  PYLITH_METHOD_END;
} // testCloneSection

// ----------------------------------------------------------------------
// Test clear().
void
pylith::topology::TestFieldMesh::testClear(void)
{ // testClear
  PYLITH_METHOD_BEGIN;

  Mesh mesh(_TestFieldMesh::cellDim);
  Field field(mesh);

  field.scale(2.0);
  field.vectorFieldType(Field::TENSOR);
  field.dimensionalizeOkay(true);
  
  field.clear();

  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0), field._metadata.scale);
  CPPUNIT_ASSERT_EQUAL(Field::OTHER, field._metadata.vectorFieldType);
  CPPUNIT_ASSERT_EQUAL(false, field._metadata.dimsOkay);

  PYLITH_METHOD_END;
} // testClear

// ----------------------------------------------------------------------
// Test allocate().
void
pylith::topology::TestFieldMesh::testAllocate(void)
{ // testAllocate
  PYLITH_METHOD_BEGIN;

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

  Field field(mesh);
  field.newSection(Field::VERTICES_FIELD, fiberDim);
  field.allocate();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d)
      fieldArray[off+d] = valuesNondim[i++];
  } // for
  fieldVisitor.clear();

  fieldVisitor.initialize(field);
  fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], fieldArray[off+d], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testAllocate

// ----------------------------------------------------------------------
// Test zero().
void
pylith::topology::TestFieldMesh::testZero(void)
{ // testZero
  PYLITH_METHOD_BEGIN;

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

  Field field(mesh);
  field.newSection(Field::VERTICES_FIELD, fiberDim);
  field.allocate();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d)
      fieldArray[off+d] = valuesNondim[i++];
  } // for
  fieldVisitor.clear();

  field.zero();

  fieldVisitor.initialize(field);
  fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, fieldArray[off+d], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testZero

// ----------------------------------------------------------------------
// Test zero().
void
pylith::topology::TestFieldMesh::testZeroAll(void)
{ // testZeroAll
  PYLITH_METHOD_BEGIN;

  const PetscInt fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };
  const PetscInt nconstraints[] = { 0, 2, 1, 3 };
  const PetscInt constraints[] = {
              // 0
    0, 2,     // 1
    2,        // 2
    0, 1, 2,  // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  // Create field and set constraint sizes
  Field field(mesh);
  field.newSection(Field::VERTICES_FIELD, fiberDim);
  PetscSection section = field.localSection();CPPUNIT_ASSERT(section);
  PetscErrorCode err = 0;
  PetscInt index = 0;
  for(PetscInt v = vStart, iV=0; v < vEnd; ++v) {
    err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]);PYLITH_CHECK_ERROR(err);
  } // for
  field.allocate();
  PetscVec vec = field.localVector();CPPUNIT_ASSERT(vec);
  for(PetscInt v = vStart, iV = 0; v < vEnd; ++v, index += nconstraints[iV++]) {
    err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[index]);PYLITH_CHECK_ERROR(err);
  } // for
  field.zero();

  VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d)
      fieldArray[off+d] = valuesNondim[i++];
  } // for
  fieldVisitor.clear();

  field.zeroAll();
  
  fieldVisitor.initialize(field);
  fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, fieldArray[off+d], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testZeroAll

// ----------------------------------------------------------------------
// Test complete().
void
pylith::topology::TestFieldMesh::testComplete(void)
{ // testComplete
  PYLITH_METHOD_BEGIN;

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

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field field(mesh);
  field.newSection(Field::VERTICES_FIELD, fiberDim);
  field.allocate();

  VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d)
      fieldArray[off+d] = valuesNondim[i++];
  } // for
  fieldVisitor.clear();

  field.complete();

  // Expect no change for this serial test
  fieldVisitor.initialize(field);
  fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], fieldArray[off+d], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testComplete

// ----------------------------------------------------------------------
// Test copy().
void
pylith::topology::TestFieldMesh::testCopy(void)
{ // testCopy
  PYLITH_METHOD_BEGIN;

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

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    VecVisitorMesh fieldVisitor(fieldSrc);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesNondim[i++];
    } // for
  } // Setup source field
    
  Field field(mesh);
  field.newSection(Field::VERTICES_FIELD, fiberDim);
  field.allocate();

  field.copy(fieldSrc);

  VecVisitorMesh fieldVisitor(field);
  const PetscScalar* fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], fieldArray[off+d], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testCopy

// ----------------------------------------------------------------------
// Test copySubfield().
void
pylith::topology::TestFieldMesh::testCopySubfield(void)
{ // testCopy
  PYLITH_METHOD_BEGIN;

  const int fiberDimA = 1;
  const char* fieldA = "one";
  const int fiberDimB = 2;
  const char* fieldB = "two";
  const int fiberDim = fiberDimA + fiberDimB;

  const PylithScalar scale = 2.0;
  const int npoints = 4;
  const PylithScalar valuesNondim[npoints*fiberDim] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };
  const PylithScalar valuesANondim[npoints*fiberDimA] = {
    1.1,
    1.2,
    1.3,
    1.4,
  };
  const PylithScalar valuesBNondim[npoints*fiberDimB] = {
    2.2, 3.3,
    2.3, 3.4,
    2.4, 3.5,
    2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field fieldSrc(mesh);
  { // Setup source field
    fieldSrc.label("solution");
    fieldSrc.subfieldAdd("one", fiberDimA, Field::SCALAR);
    fieldSrc.subfieldAdd("two", fiberDimB, Field::VECTOR);
    fieldSrc.subfieldsSetup();
    fieldSrc.newSection(Field::VERTICES_FIELD, fiberDim);
    fieldSrc.subfieldSetDof("one", Field::VERTICES_FIELD, fiberDimA);
    fieldSrc.subfieldSetDof("two", Field::VERTICES_FIELD, fiberDimB);
    fieldSrc.allocate();

    VecVisitorMesh fieldVisitor(fieldSrc);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for (PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
      for (PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesNondim[i++];
    } // for
  } // Setup source field
    
  { // Test with preallocated field
    Field field(mesh);
    field.newSection(Field::VERTICES_FIELD, fiberDimB);
    field.allocate();
    field.copySubfield(fieldSrc, "two");

    VecVisitorMesh fieldVisitor(field);
    const PetscScalar* fieldArray = fieldVisitor.localArray();
    const PylithScalar tolerance = 1.0e-6;
    for (PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(fiberDimB, fieldVisitor.sectionDof(v));
      for (PetscInt d = 0; d < fiberDimB; ++d) {
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesBNondim[i++], fieldArray[off+d], tolerance);
      } // for
    } // for
  } // Test with preallocated field

  { // Test with unallocated field
    Field field(mesh);
    field.copySubfield(fieldSrc, "one");

    VecVisitorMesh fieldVisitor(field);
    const PetscScalar* fieldArray = fieldVisitor.localArray();
    const PylithScalar tolerance = 1.0e-6;
    for (PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(fiberDimA, fieldVisitor.sectionDof(v));
      for (PetscInt d = 0; d < fiberDimA; ++d) {
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesANondim[i++], fieldArray[off+d], tolerance);
      } // for
    } // for
  } // Test with unallocated field

  PYLITH_METHOD_END;
} // testCopySubfield

// ----------------------------------------------------------------------
// Test operator+=().
void
pylith::topology::TestFieldMesh::testOperatorAdd(void)
{ // testOperateAdd
  PYLITH_METHOD_BEGIN;

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

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    VecVisitorMesh fieldVisitor(fieldSrc);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesA[i++];
    } // for
  } // Setup source field
    

  Field field(mesh);
  { // Setup destination field
    field.newSection(Field::VERTICES_FIELD, fiberDim);
    field.allocate();
    VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesB[i++];
    } // for
  } // Setup destination field

  field += fieldSrc;

  VecVisitorMesh fieldVisitor(field);
  const PetscScalar* fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d, ++i) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesA[i] + valuesB[i], fieldArray[off+d], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testOperateAdd

// ----------------------------------------------------------------------
// Test dimensionalize().
void
pylith::topology::TestFieldMesh::testDimensionalize(void)
{ // testDimensionalize
  PYLITH_METHOD_BEGIN;

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

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field field(mesh);
  { // setup field
    field.newSection(Field::VERTICES_FIELD, fiberDim);
    field.allocate();

    VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesNondim[i++];
    } // for
  } // setup field

  field.scale(scale);
  field.dimensionalizeOkay(true);
  field.dimensionalize();

  VecVisitorMesh fieldVisitor(field);
  const PetscScalar* fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d, ++i) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i]*scale, fieldArray[off+d], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testDimensionalize

// ----------------------------------------------------------------------
// Test view().
void
pylith::topology::TestFieldMesh::testView(void)
{ // testView
  PYLITH_METHOD_BEGIN;

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

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field field(mesh);
  field.newSection(Field::VERTICES_FIELD, fiberDim);
  field.allocate();
  VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d)
      fieldArray[off+d] = valuesNondim[i++];
  } // for

  field.view("Testing view");

  PYLITH_METHOD_END;
} // testView

// ----------------------------------------------------------------------
// Test createScatter().
void
pylith::topology::TestFieldMesh::testCreateScatter(void)
{ // testCreateScatter
  PYLITH_METHOD_BEGIN;

  const int fiberDim = 3;

  Mesh mesh;
  _buildMesh(&mesh);

  Field field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(Field::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatter(mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  const Field::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.dm);
  CPPUNIT_ASSERT(sinfo.vector);

  const int sizeE = (vEnd-vStart) * fiberDim;
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
  const Field::ScatterInfo& sinfoB = field._getScatter("B");
  CPPUNIT_ASSERT(sinfoB.dm);
  CPPUNIT_ASSERT(sinfoB.vector);

  Field field2(mesh);
  field2.cloneSection(field);
  CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

  const Field::ScatterInfo& sinfo2 = field2._getScatter("");
  CPPUNIT_ASSERT(sinfo2.dm);
  CPPUNIT_ASSERT(sinfo2.vector);

  const Field::ScatterInfo& sinfo2B = field2._getScatter("B");
  CPPUNIT_ASSERT(sinfo2B.dm);
  CPPUNIT_ASSERT(sinfo2B.vector);


  PYLITH_METHOD_END;
} // testCreateScatter

// ----------------------------------------------------------------------
// Test createScatterWithBC().
void
pylith::topology::TestFieldMesh::testCreateScatterWithBC(void)
{ // testCreateScatterWithBC
  PYLITH_METHOD_BEGIN;

  const int fiberDim = 3;

  Mesh mesh;
  _buildMesh(&mesh);

  Field field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(Field::VERTICES_FIELD, fiberDim);
  field.allocate();

  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatterWithBC(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  const Field::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.dm);
  CPPUNIT_ASSERT(sinfo.vector);

  const int sizeE = (vEnd-vStart) * fiberDim;
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
  const Field::ScatterInfo& sinfoB = field._getScatter("B");
  CPPUNIT_ASSERT(sinfoB.dm);
  CPPUNIT_ASSERT(sinfoB.vector);

  Field field2(mesh);
  field2.cloneSection(field);
  CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

  const Field::ScatterInfo& sinfo2 = field2._getScatter("");
  CPPUNIT_ASSERT(sinfo2.dm);
  CPPUNIT_ASSERT(sinfo2.vector);

  const Field::ScatterInfo& sinfo2B = field2._getScatter("B");
  CPPUNIT_ASSERT(sinfo2B.dm);
  CPPUNIT_ASSERT(sinfo2B.vector);


  PYLITH_METHOD_END;
} // testCreateScatterWithBC

// ----------------------------------------------------------------------
// Test vector().
void
pylith::topology::TestFieldMesh::testVector(void)
{ // testVector
  PYLITH_METHOD_BEGIN;

  const int fiberDim = 3;

  Mesh mesh;
  _buildMesh(&mesh);

  Field field(mesh);
  field.newSection(Field::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatter(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.dm);
  CPPUNIT_ASSERT(sinfo.vector);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  const PetscVec vec = field.vector();
  CPPUNIT_ASSERT_EQUAL(sinfo.vector, vec);
  int size = 0;
  PetscErrorCode err = VecGetSize(vec, &size);PYLITH_CHECK_ERROR(err);
  const int sizeE = (vEnd-vStart) * fiberDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  PYLITH_METHOD_END;
} // testVector

// ----------------------------------------------------------------------
// Test scatterLocalToGlobal().
void
pylith::topology::TestFieldMesh::testScatterLocalToGlobal(void)
{ // testScatterLocalToGlobal
  PYLITH_METHOD_BEGIN;

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

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field field(mesh);
  { // setup field
    field.newSection(Field::VERTICES_FIELD, fiberDim);
    field.allocate();
    VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesE[i++];
    } // for
  } // setup field

  field.createScatter(mesh, context);
  field.scatterLocalToGlobal(context);

  const PetscVec vec = field.vector(context);CPPUNIT_ASSERT(vec);
  PetscInt size = 0;
  PetscErrorCode err = VecGetSize(vec, &size);PYLITH_CHECK_ERROR(err);
  PetscScalar* valuesVec = NULL;
  err = VecGetArray(vec, &valuesVec);PYLITH_CHECK_ERROR(err);

  const PylithScalar tolerance = 1.0e-06;
  const int sizeE = (vEnd-vStart) * fiberDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);
  for (int i=0; i < sizeE; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i], valuesVec[i], tolerance);
  err = VecRestoreArray(vec, &valuesVec);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // testScatterLocalToGlobal

// ----------------------------------------------------------------------
// Test scatterGlobalToLocal().
void
pylith::topology::TestFieldMesh::testScatterGlobalToLocal(void)
{ // testScatterGlobalToLocal
  PYLITH_METHOD_BEGIN;

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

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field field(mesh);
  field.newSection(Field::VERTICES_FIELD, fiberDim);
  field.allocate();
  field.createScatter(mesh, context);

  const PetscVec vec = field.vector(context);CPPUNIT_ASSERT(vec);
  PetscInt size = 0;
  PetscErrorCode err = VecGetSize(vec, &size);PYLITH_CHECK_ERROR(err);
  PetscScalar* valuesVec = NULL;
  err = VecGetArray(vec, &valuesVec);PYLITH_CHECK_ERROR(err);
  const int sizeE = (vEnd-vStart) * fiberDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);
  for (int i=0; i < sizeE; ++i)
    valuesVec[i] = valuesE[i];
  err = VecRestoreArray(vec, &valuesVec);PYLITH_CHECK_ERROR(err);

  field.scatterGlobalToLocal(context);

  const PylithScalar tolerance = 1.0e-06;
  VecVisitorMesh fieldVisitor(field);
  const PetscScalar* fieldArray = fieldVisitor.localArray();
  for (PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i++], fieldArray[off+iDim], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testScatterGlobalToLocal

// ----------------------------------------------------------------------
// Test splitDefault().
void
pylith::topology::TestFieldMesh::testSplitDefault(void)
{ // testSplitDefault
  PYLITH_METHOD_BEGIN;

  const PetscInt spaceDim = _TestFieldMesh::cellDim;
  const PetscInt numFields = spaceDim;
  const PetscInt nconstraints[4] = { 1, 2, 0, 1 };
  const PetscInt constraints[4] = {
    1,     // 0
    0, 1,  // 1
           // 2
    0,     // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  PetscErrorCode err = 0;

  // Create field with section to use to create new field
  Field fieldSrc(mesh);
  { // Setup source field
    for(PetscInt f = 0; f < numFields; ++f) {
      std::ostringstream msg;
      msg << "Field "<<f;
      fieldSrc.subfieldAdd(msg.str().c_str(), 1, Field::SCALAR);
    } // for
    fieldSrc.subfieldsSetup();
    fieldSrc.newSection(Field::VERTICES_FIELD, spaceDim);
    for(PetscInt f = 0; f < spaceDim; ++f) {
      std::ostringstream msg;
      msg << "Field "<<f;
      fieldSrc.subfieldSetDof(msg.str().c_str(), Field::VERTICES_FIELD, 1);
    } // for
    PetscSection section = fieldSrc.localSection();CPPUNIT_ASSERT(section);
    PetscInt iV = 0, iC = 0;
    for(PetscInt v = vStart; v < vEnd; ++v) {
      const int nconstraintsVertex = nconstraints[iV];
      err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]);PYLITH_CHECK_ERROR(err);
      for(PetscInt iConstraint = 0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const int field = constraints[iC++];
        err = PetscSectionAddFieldConstraintDof(section, v, field, 1);PYLITH_CHECK_ERROR(err);
      } // for
    } // for
    fieldSrc.allocate();

    iV = 0; iC = 0;
    PetscInt zero = 0;
    for(PetscInt v = vStart; v < vEnd; ++v) {
      const int nconstraintsVertex = nconstraints[iV++];
      if (nconstraintsVertex > 0) {
        err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[iC]);PYLITH_CHECK_ERROR(err);
      } // if
      for(PetscInt iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const PetscInt field = constraints[iC++];
        err = PetscSectionSetFieldConstraintIndices(section, v, field, &zero);PYLITH_CHECK_ERROR(err);
      } // for
    } // for
  } // Setup source field

  PetscSection section = fieldSrc.localSection();CPPUNIT_ASSERT(section);
  PetscInt numSectionFields;
  err = PetscSectionGetNumFields(section, &numSectionFields);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(numFields, numSectionFields);

  for(PetscInt f = 0; f < numFields; ++f) {
    PetscInt iV = 0, iC = 0;

    for(PetscInt v = vStart; v < vEnd; ++v, ++iV) {
      PetscInt fdof, fcdof;
      PetscBool isConstrained = PETSC_FALSE;
      const PetscInt nconstraintsVertex = nconstraints[iV];

      err = PetscSectionGetFieldDof(section, v, f, &fdof);PYLITH_CHECK_ERROR(err);
      err = PetscSectionGetFieldConstraintDof(section, v, f, &fcdof);PYLITH_CHECK_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(1, fdof);
      for(PetscInt iConstraint = 0; iConstraint < nconstraintsVertex; ++iConstraint)
        if (constraints[iC++] == f) isConstrained = PETSC_TRUE;
      const PetscInt constraintDimE = (!isConstrained) ? 0 : 1;
      CPPUNIT_ASSERT_EQUAL(constraintDimE, fcdof);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testSplitDefault

// ----------------------------------------------------------------------
// Test cloneSection() with split field.
void
pylith::topology::TestFieldMesh::testCloneSectionSplit(void)
{ // testCloneSectionSplit
  PYLITH_METHOD_BEGIN;

  const PetscInt spaceDim = _TestFieldMesh::cellDim;
  const PetscInt numFields = spaceDim;
  const PetscInt nconstraints[4] = { 1, 2, 0, 1 };
  const PetscInt constraints[4] = {
    1,     // 0
    0, 1,  // 1
           // 2
    0,     // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  PetscErrorCode err = 0;

  // Create field with atlas to use to create new field
  Field fieldSrc(mesh);
  { // Setup source field
    for(PetscInt f = 0; f < numFields; ++f) {
      std::ostringstream msg;
      msg << "Field "<<f;
      fieldSrc.subfieldAdd(msg.str().c_str(), 1, Field::SCALAR);
    } // for
    fieldSrc.subfieldsSetup();
    fieldSrc.newSection(Field::VERTICES_FIELD, spaceDim);
    for(PetscInt f = 0; f < spaceDim; ++f) {
      std::ostringstream msg;
      msg << "Field "<<f;
      fieldSrc.subfieldSetDof(msg.str().c_str(), Field::VERTICES_FIELD, 1);
    } // for
    PetscSection section = fieldSrc.localSection();
    CPPUNIT_ASSERT(section);
    PetscInt iV = 0, iC = 0;
    for(PetscInt v = vStart; v < vEnd; ++v) {
      const int nconstraintsVertex = nconstraints[iV];
      err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]);PYLITH_CHECK_ERROR(err);
      for(PetscInt iConstraint = 0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const int field = constraints[iC++];
        err = PetscSectionAddFieldConstraintDof(section, v, field, 1);PYLITH_CHECK_ERROR(err);
      } // for
    } // for
    fieldSrc.allocate();

    iV = 0; iC = 0;
    PetscInt zero = 0;
    for(PetscInt v = vStart; v < vEnd; ++v) {
      const int nconstraintsVertex = nconstraints[iV++];
      if (nconstraintsVertex > 0) {
        err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[iC]);PYLITH_CHECK_ERROR(err);
      } // if
      for(PetscInt iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const PetscInt field = constraints[iC++];
        err = PetscSectionSetFieldConstraintIndices(section, v, field, &zero);PYLITH_CHECK_ERROR(err);
      } // for
    } // for
  } // Setup source field

  Field field(mesh);
  field.cloneSection(fieldSrc);

  PetscSection section = field.localSection();CPPUNIT_ASSERT(section);
  PetscInt numSectionFields;
  err = PetscSectionGetNumFields(section, &numSectionFields);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(numFields, numSectionFields);

  for(PetscInt f = 0; f < numFields; ++f) {
    PetscInt iV = 0, iC = 0;

    for(PetscInt v = vStart; v < vEnd; ++v, ++iV) {
      PetscInt fdof, fcdof;
      PetscBool isConstrained = PETSC_FALSE;
      const PetscInt nconstraintsVertex = nconstraints[iV];

      err = PetscSectionGetFieldDof(section, v, f, &fdof);PYLITH_CHECK_ERROR(err);
      err = PetscSectionGetFieldConstraintDof(section, v, f, &fcdof);PYLITH_CHECK_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(1, fdof);
      for(PetscInt iConstraint = 0; iConstraint < nconstraintsVertex; ++iConstraint)
        if (constraints[iC++] == f) isConstrained = PETSC_TRUE;
      const PetscInt constraintDimE = (!isConstrained) ? 0 : 1;
      CPPUNIT_ASSERT_EQUAL(constraintDimE, fcdof);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testCloneSectionSplit

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_buildMesh(Mesh* mesh)
{ // _buildMesh
  PYLITH_METHOD_BEGIN;

  assert(mesh);

  const int cellDim = _TestFieldMesh::cellDim;
  const int ncells = _TestFieldMesh::ncells;
  const int* cells = _TestFieldMesh::cells;
  const int nvertices = _TestFieldMesh::nvertices;
  const int ncorners = _TestFieldMesh::ncorners;
  const int spaceDim = _TestFieldMesh::cellDim;
  const PylithScalar* coordinates = _TestFieldMesh::coordinates;
  const bool interpolate = false;

  PetscErrorCode err = 0;

  MeshOps::createDMMesh(mesh, _TestFieldMesh::cellDim);
  PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  
  err = DMPlexSetChart(dmMesh, 0, ncells+nvertices);PYLITH_CHECK_ERROR(err);
  for(PetscInt c = 0; c < ncells; ++c) {
    err = DMPlexSetConeSize(dmMesh, c, ncorners);PYLITH_CHECK_ERROR(err);
  } // for
  err = DMSetUp(dmMesh);PYLITH_CHECK_ERROR(err);
  PetscInt *cone = new PetscInt[ncorners];
  for(PetscInt c = 0; c < ncells; ++c) {
    for(PetscInt v = 0; v < ncorners; ++v) {
      cone[v] = cells[c*ncorners+v]+ncells;
    } // for
    err = DMPlexSetCone(dmMesh, c, cone);PYLITH_CHECK_ERROR(err);
  } // for
  delete[] cone; cone = 0;
  err = DMPlexSymmetrize(dmMesh);PYLITH_CHECK_ERROR(err);
  err = DMPlexStratify(dmMesh);PYLITH_CHECK_ERROR(err);
  PetscSection coordSection = NULL;
  PetscVec coordVec = NULL;
  PetscScalar *coords = NULL;
  PetscInt coordSize;

  err = DMGetCoordinateSection(dmMesh, &coordSection);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetNumFields(coordSection, 1);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetFieldComponents(coordSection, 0, spaceDim);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetChart(coordSection, ncells, ncells+nvertices);PYLITH_CHECK_ERROR(err);
  for(PetscInt v = ncells; v < ncells+nvertices; ++v) {
    err = PetscSectionSetDof(coordSection, v, spaceDim);PYLITH_CHECK_ERROR(err);
  } // for
  err = PetscSectionSetUp(coordSection);PYLITH_CHECK_ERROR(err);
  err = PetscSectionGetStorageSize(coordSection, &coordSize);PYLITH_CHECK_ERROR(err);
  err = VecCreate(mesh->comm(), &coordVec);PYLITH_CHECK_ERROR(err);
  err = VecSetSizes(coordVec, coordSize, PETSC_DETERMINE);PYLITH_CHECK_ERROR(err);
  err = VecSetFromOptions(coordVec);PYLITH_CHECK_ERROR(err);
  err = VecGetArray(coordVec, &coords);PYLITH_CHECK_ERROR(err);
  for(PetscInt v = 0; v < nvertices; ++v) {
    PetscInt off;
    err = PetscSectionGetOffset(coordSection, v+ncells, &off);PYLITH_CHECK_ERROR(err);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      coords[off+d] = coordinates[v*spaceDim+d];
    } // for
  } // for
  err = VecRestoreArray(coordVec, &coords);PYLITH_CHECK_ERROR(err);
  err = DMSetCoordinatesLocal(dmMesh, coordVec);PYLITH_CHECK_ERROR(err);
  err = VecDestroy(&coordVec);PYLITH_CHECK_ERROR(err);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);

  PYLITH_METHOD_END;
} // _buildMesh


// End of file 
