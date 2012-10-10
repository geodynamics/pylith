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
  PetscSection section = field.petscSection();
  CPPUNIT_ASSERT(!section);
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());

  field.newSection();
  PetscSection section = field.petscSection();
  CPPUNIT_ASSERT(section);
} // testNewSection

// ----------------------------------------------------------------------
// Test newSection(points).
void
pylith::topology::TestFieldMesh::testNewSectionPoints(void)
{ // testNewSectionPoints
  const int fiberDim = 2;

  Mesh mesh;
  _buildMesh(&mesh);
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);

  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof;
    err = PetscSectionGetDof(section, v, &dof);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDim, dof);
  }
  const char *name;
  err = PetscObjectGetName((PetscObject) vec, &name);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));
} // testNewSectionPoints

// ----------------------------------------------------------------------
// Test newSection(int_array).
void
pylith::topology::TestFieldMesh::testNewSectionPointsArray(void)
{ // testNewSectionPointsArray
  const int fiberDim = 2;

  Mesh mesh;
  _buildMesh(&mesh);
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
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
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);

  // Points in array should have a fiber dimension of fiberDim.
  for(int i = 0; i < pointsIn.size(); ++i) {
    PetscInt dof;
    err = PetscSectionGetDof(section, pointsIn[i], &dof);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDim, dof);
  }

  // Points not int array should have a fiber dimension of zero.
  PetscInt pStart, pEnd;
  err = PetscSectionGetChart(section, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  for(int i = 0; i < pointsOut.size(); ++i) {
    if (pointsOut[i] >= pStart && pointsOut[i] < pEnd) {
      PetscInt dof;
      err = PetscSectionGetDof(section, pointsOut[i], &dof);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(0, dof);
    }
  }
  const char *name;
  err = PetscObjectGetName((PetscObject) vec, &name);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));
} // testNewSectionPointsArray

// ----------------------------------------------------------------------
// Test newSection(domain).
void
pylith::topology::TestFieldMesh::testNewSectionDomain(void)
{ // testNewSectionDomain
  const int fiberDim = 2;

  Mesh mesh;
  _buildMesh(&mesh);
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof;
    err = PetscSectionGetDof(section, v, &dof);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDim, dof);
  }
  const char *name;
  err = PetscObjectGetName((PetscObject) vec, &name);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));
} // testNewSectionDomain

// ----------------------------------------------------------------------
// Test newSection(field).
void
pylith::topology::TestFieldMesh::testNewSectionField(void)
{ // testNewSectionField
  const int fiberDim = 3;
    
  Mesh mesh;
  _buildMesh(&mesh);
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  fieldSrc.allocate();

  const int fiberDim2 = 5;
  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(fieldSrc, fiberDim2);
  field.allocate();
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof;
    err = PetscSectionGetDof(section, v, &dof);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDim2, dof);
  }
  const char *name;
  err = PetscObjectGetName((PetscObject) vec, &name);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));
} // testNewSectionField

// ----------------------------------------------------------------------
// Test cloneSection().
void
pylith::topology::TestFieldMesh::testCloneSection(void)
{ // testCloneSection
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    PetscSection section = fieldSrc.petscSection();
    CPPUNIT_ASSERT(section);
    int iV=0;
    for(PetscInt v = vStart; v < vEnd; ++v) {
      err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]);CHECK_PETSC_ERROR(err);
    }
    fieldSrc.allocate();

    int index = 0;
    iV = 0;
    for(PetscInt v = vStart; v < vEnd; ++v, index += nconstraints[iV++]) {
      err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[index]);CHECK_PETSC_ERROR(err);
    }
    fieldSrc.zero();
    fieldSrc.createScatter(mesh);
    fieldSrc.createScatter(mesh, "A");
  } // Setup source field

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.cloneSection(fieldSrc);
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
  int iV = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, cdof;
    err = PetscSectionGetDof(section, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintDof(section, v, &cdof);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDim, dof);
    CPPUNIT_ASSERT_EQUAL(nconstraints[iV++], cdof);
  } // for

  // Verify vector scatters were also copied.
  CPPUNIT_ASSERT_EQUAL(fieldSrc._scatters[""].dm,  field._scatters[""].dm);
  CPPUNIT_ASSERT_EQUAL(fieldSrc._scatters["A"].dm, field._scatters["A"].dm);
  const char *name;
  err = PetscObjectGetName((PetscObject) vec, &name);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);

  PetscScalar *array;
  PetscInt     i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d)
      array[off+d] = valuesNondim[i++];
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);

  const PylithScalar tolerance = 1.0e-6;
  i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], array[off+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);

  PetscScalar *array;
  PetscInt     i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d)
      array[off+d] = valuesNondim[i++];
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);

  field.zero();

  const PylithScalar tolerance = 1.0e-6;
  i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, array[off+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
} // testZero

// ----------------------------------------------------------------------
// Test zero().
void
pylith::topology::TestFieldMesh::testZeroAll(void)
{ // testZeroAll
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  // Create field and set constraint sizes
  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  PetscSection section = field.petscSection();
  PetscInt     iV      = 0;
  PetscInt     index   = 0;
  CPPUNIT_ASSERT(section);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]);CHECK_PETSC_ERROR(err);
  }
  field.allocate();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(vec);
  iV = 0;
  for(PetscInt v = vStart; v < vEnd; ++v, index += nconstraints[iV++]) {
    err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[index]);CHECK_PETSC_ERROR(err);
  }
  field.zero();

  PetscScalar *array;
  PetscInt     i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d)
      array[off+d] = valuesNondim[i++];
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
  
  field.zeroAll();
  
  const PylithScalar tolerance = 1.0e-6;
  i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, array[off+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);

  PetscScalar *array;
  PetscInt     i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d)
      array[off+d] = valuesNondim[i++];
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);

  field.complete();

  // Expect no change for this serial test
  const PylithScalar tolerance = 1.0e-6;
  i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], array[off+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    PetscSection section = fieldSrc.petscSection();
    Vec          vec     = fieldSrc.localVector();
    CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
    
    PetscScalar *array;
    PetscInt     i = 0;
    err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v) {
      PetscInt off;

      err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < fiberDim; ++d)
        array[off+d] = valuesNondim[i++];
    } // for
    err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
  } // Setup source field

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);

  field.copy(fieldSrc);

  const PylithScalar tolerance = 1.0e-6;
  PetscScalar *array;
  PetscInt i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], array[off+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    PetscSection section = fieldSrc.petscSection();
    Vec          vec     = fieldSrc.localVector();
    CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
    
    PetscScalar *array;
    PetscInt     i = 0;
    err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v) {
      PetscInt off;

      err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < fiberDim; ++d)
        array[off+d] = valuesA[i++];
    } // for
    err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
  } // Setup source field

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
  { // Setup destination field

    PetscScalar *array;
    PetscInt     i = 0;
    err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v) {
      PetscInt off;

      err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < fiberDim; ++d)
        array[off+d] = valuesB[i++];
    } // for
    err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
  } // Setup destination field

  field += fieldSrc;

  const PylithScalar tolerance = 1.0e-6;
  PetscScalar *array;
  PetscInt i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesA[i] + valuesB[i], array[off+d], tolerance);
      ++i;
    } // for
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);

  PetscScalar *array;
  PetscInt     i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d)
      array[off+d] = valuesNondim[i++];
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);

  field.scale(scale);
  field.addDimensionOkay(true);
  field.dimensionalize();

  const PylithScalar tolerance = 1.0e-6;
  i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++]*scale, array[off+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);

  PetscScalar *array;
  PetscInt     i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d)
      array[off+d] = valuesNondim[i++];
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);

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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  const int sizeE = (vEnd-vStart) * fiberDim;

  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatter(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field<Mesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.dm);
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
  CPPUNIT_ASSERT(sinfoB.dm);
  CPPUNIT_ASSERT(sinfoB.vector);

  Field<Mesh> field2(mesh);
  field2.cloneSection(field);
  CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

  const Field<Mesh>::ScatterInfo& sinfo2 = field2._getScatter("");
  CPPUNIT_ASSERT(sinfo2.dm);
  CPPUNIT_ASSERT(sinfo2.vector);

  const Field<Mesh>::ScatterInfo& sinfo2B = field2._getScatter("B");
  CPPUNIT_ASSERT(sinfo2B.dm);
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  const int sizeE = (vEnd-vStart) * fiberDim;

  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatterWithBC(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field<Mesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.dm);
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
  CPPUNIT_ASSERT(sinfoB.dm);
  CPPUNIT_ASSERT(sinfoB.vector);

  Field<Mesh> field2(mesh);
  field2.cloneSection(field);
  CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

  const Field<Mesh>::ScatterInfo& sinfo2 = field2._getScatter("");
  CPPUNIT_ASSERT(sinfo2.dm);
  CPPUNIT_ASSERT(sinfo2.vector);

  const Field<Mesh>::ScatterInfo& sinfo2B = field2._getScatter("B");
  CPPUNIT_ASSERT(sinfo2B.dm);
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatter(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field<Mesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.dm);
  CPPUNIT_ASSERT(sinfo.vector);

  const PetscVec vec = field.vector();
  CPPUNIT_ASSERT_EQUAL(sinfo.vector, vec);
  int size = 0;
  err = VecGetSize(vec, &size);CHECK_PETSC_ERROR(err);
  const int sizeE = (vEnd-vStart) * fiberDim;
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  PetscSection section = field.petscSection();
  Vec          secvec  = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(secvec);

  PetscScalar *array;
  PetscInt     i = 0;
  err = VecGetArray(secvec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d)
      array[off+d] = valuesE[i++];
  } // for
  err = VecRestoreArray(secvec, &array);CHECK_PETSC_ERROR(err);

  field.createScatter(mesh, context);
  field.scatterSectionToVector(context);
  const PetscVec vec = field.vector(context);
  CPPUNIT_ASSERT(0 != vec);
  PetscInt size = 0;
  err = VecGetSize(vec, &size);CHECK_PETSC_ERROR(err);
  PylithScalar* valuesVec = 0;
  err = VecGetArray(vec, &valuesVec);CHECK_PETSC_ERROR(err);

  const PylithScalar tolerance = 1.0e-06;
  const int sizeE = (vEnd-vStart) * fiberDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);
  for (int i=0; i < sizeE; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i], valuesVec[i], tolerance);
  err = VecRestoreArray(vec, &valuesVec);CHECK_PETSC_ERROR(err);
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  PetscSection section = field.petscSection();
  CPPUNIT_ASSERT(section);

  field.createScatter(mesh, context);
  const PetscVec vec = field.vector(context);
  CPPUNIT_ASSERT(0 != vec);
  PetscInt size = 0;
  err = VecGetSize(vec, &size);CHECK_PETSC_ERROR(err);
  PylithScalar* valuesVec = 0;
  err = VecGetArray(vec, &valuesVec);CHECK_PETSC_ERROR(err);

  const PylithScalar tolerance = 1.0e-06;
  const int sizeE = (vEnd-vStart) * fiberDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);
  for (int i=0; i < sizeE; ++i)
    valuesVec[i] = valuesE[i];
  err = VecRestoreArray(vec, &valuesVec);CHECK_PETSC_ERROR(err);

  field.createScatter(mesh, context);
  field.scatterVectorToSection(context);

  PetscScalar *array;
  PetscInt     i = 0;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < fiberDim; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i++], array[off+d], tolerance);
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  CPPUNIT_ASSERT(0);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, spaceDim);
    fieldSrc.splitDefault();
    const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
    CPPUNIT_ASSERT(!section.isNull());
    int iV=0;
    int iC=0;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iV) {
      const int nconstraintsVertex = nconstraints[iV];
      section->addConstraintDimension(v, nconstraintsVertex);
      for (int iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const int fibration = constraints[iC++];
        section->addConstraintDimension(v, 1, fibration);
      } // for
    } // for
    fieldSrc.allocate();

    iC = 0;
    iV = 0;
    int zero = 0;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iV) {
      const int nconstraintsVertex = nconstraints[iV];
      if (nconstraintsVertex > 0)
        section->setConstraintDof(v, &constraints[iC]);
      for (int iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const int fibration = constraints[iC++];
        section->setConstraintDof(v, &zero, fibration);
      } // for
    } // for
  } // Setup source field

  const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
  CPPUNIT_ASSERT(!section.isNull());
  CPPUNIT_ASSERT_EQUAL(numFibrations, section->getNumSpaces());

  for (int fibration=0; fibration < spaceDim; ++fibration) {
    const ALE::Obj<Mesh::RealSection>& sectionSplit = section->getFibration(fibration);
    CPPUNIT_ASSERT(!sectionSplit.isNull());
    int iV = 0;
    int iC = 0;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iV) {
      CPPUNIT_ASSERT_EQUAL(1, section->getFiberDimension(v, fibration));
      bool isConstrained = false;
      const int nconstraintsVertex = nconstraints[iV];
      for (int iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint)
        if (constraints[iC++] == fibration)
          isConstrained = true;
      const int constraintDimE = (!isConstrained) ? 0 : 1;
      CPPUNIT_ASSERT_EQUAL(constraintDimE,
                           section->getConstraintDimension(v, fibration));
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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  CPPUNIT_ASSERT(0);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, spaceDim);
    fieldSrc.splitDefault();
    const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
    CPPUNIT_ASSERT(!section.isNull());
    int iV=0;
    int iC=0;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iV) {
      const int nconstraintsVertex = nconstraints[iV];
      section->addConstraintDimension(v, nconstraintsVertex);
      for (int iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const int fibration = constraints[iC++];
        section->addConstraintDimension(v, 1, fibration);
      } // for
    } // for
    fieldSrc.allocate();

    iC = 0;
    iV = 0;
    int zero = 0;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iV) {
      const int nconstraintsVertex = nconstraints[iV];
      if (nconstraintsVertex > 0)
        section->setConstraintDof(v, &constraints[iC]);
      for (int iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const int fibration = constraints[iC++];
        section->setConstraintDof(v, &zero, fibration);
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
    int iV = 0;
    int iC = 0;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iV) {
      CPPUNIT_ASSERT_EQUAL(1, section->getFiberDimension(v, fibration));
      bool isConstrained = false;
      const int nconstraintsVertex = nconstraints[iV];
      for (int iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint)
        if (constraints[iC++] == fibration)
          isConstrained = true;
      const int constraintDimE = (!isConstrained) ? 0 : 1;
      CPPUNIT_ASSERT_EQUAL(constraintDimE,
                           section->getConstraintDimension(v, fibration));
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

  mesh->createDMMesh(_TestFieldMesh::cellDim);
  DM dmMesh = mesh->dmMesh();
  PetscErrorCode err;
  
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

  err = DMComplexSetChart(dmMesh, 0, ncells+nvertices);CHECK_PETSC_ERROR(err);
  for(PetscInt c = 0; c < ncells; ++c) {
    err = DMComplexSetConeSize(dmMesh, c, ncorners);CHECK_PETSC_ERROR(err);
  }
  err = DMSetUp(dmMesh);CHECK_PETSC_ERROR(err);
  PetscInt *cone = new PetscInt[ncorners];
  for(PetscInt c = 0; c < ncells; ++c) {
    for(PetscInt v = 0; v < ncorners; ++v) {
      cone[v] = cells[c*ncorners+v]+ncells;
    }
    err = DMComplexSetCone(dmMesh, c, cone);CHECK_PETSC_ERROR(err);
  } // for
  delete [] cone; cone = 0;
  err = DMComplexSymmetrize(dmMesh);CHECK_PETSC_ERROR(err);
  err = DMComplexStratify(dmMesh);CHECK_PETSC_ERROR(err);
  PetscSection coordSection;
  Vec          coordVec;
  PetscScalar *coords;
  PetscInt     coordSize;

  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(coordSection, ncells, ncells+nvertices);CHECK_PETSC_ERROR(err);
  for(PetscInt v = ncells; v < ncells+nvertices; ++v) {
    err = PetscSectionSetDof(coordSection, v, spaceDim);CHECK_PETSC_ERROR(err);
  }
  err = PetscSectionSetUp(coordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetStorageSize(coordSection, &coordSize);CHECK_PETSC_ERROR(err);
  err = VecCreate(sieveMesh->comm(), &coordVec);CHECK_PETSC_ERROR(err);
  err = VecSetSizes(coordVec, coordSize, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(coordVec);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  for(PetscInt v = 0; v < nvertices; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(coordSection, v+ncells, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      coords[off+d] = coordinates[v*spaceDim+d];
    }
  }
  err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  err = DMSetCoordinatesLocal(dmMesh, coordVec);CHECK_PETSC_ERROR(err);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);
} // _buildMesh


// End of file 
