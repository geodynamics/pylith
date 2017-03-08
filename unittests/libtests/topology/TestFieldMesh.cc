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
// Copyright (c) 2010-2016 University of California, Davis
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

#include "pylith/meshio/MeshBuilder.hh" // Uses MeshBuilder
#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::topology::TestFieldMesh::setUp(void)
{ // setUp
    PYLITH_METHOD_BEGIN;

    _data = new TestFieldMesh_Data; CPPUNIT_ASSERT(_data);

    PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::topology::TestFieldMesh::tearDown(void)
{ // tearDown
    PYLITH_METHOD_BEGIN;

    delete _data; _data = NULL;

    PYLITH_METHOD_END;
} // tearDown

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
// Test mesh().
void
pylith::topology::TestFieldMesh::testMesh(void)
{ // testMesh
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh mesh;
    _initializeMesh(&mesh);
    Field field(mesh);

    CPPUNIT_ASSERT_EQUAL(_data->cellDim, field.mesh().dimension());

    PYLITH_METHOD_END;
} // testMesh

// ----------------------------------------------------------------------
// Test label(), vectorFieldType(), scale(), addDimensionOkay(), spaceDim().
void
pylith::topology::TestFieldMesh::testGeneralAccessors(void)
{ // testGeneralAccessors
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh mesh;
    _initializeMesh(&mesh);
    Field field(mesh);

    // Test label()
    const std::string label = "velocity";
    field.label(label.c_str());
    CPPUNIT_ASSERT_EQUAL(label, std::string(field.label()));

    // Test vectorFieldType()
    field.vectorFieldType(FieldBase::SCALAR);
    CPPUNIT_ASSERT_EQUAL(FieldBase::SCALAR, field._metadata.vectorFieldType);

    // Test scale()
    const PylithScalar scale = 2.0;
    field.scale(scale);
    const PylithScalar tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(scale, field._metadata.scale, tolerance);

    // Test addDimensionOkay()
    CPPUNIT_ASSERT_EQUAL(false, field._metadata.dimsOkay);
    field.dimensionalizeOkay(true);
    CPPUNIT_ASSERT_EQUAL(true, field._metadata.dimsOkay);

    // Test spaceDim()
    CPPUNIT_ASSERT_EQUAL(_data->cellDim, field.spaceDim());

    PYLITH_METHOD_END;
} // testGeneralAccessors

// ----------------------------------------------------------------------
// Test chartSize(), sectionSize(), localSection(), globalSection().
void
pylith::topology::TestFieldMesh::testSectionAccessors(void)
{ // testSectionAccessors
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh mesh;
    _initializeMesh(&mesh);

    const int fiberDim = 2;
    const std::string& label = "field A";
    Field field(mesh);
    field.label(label.c_str());

    // Tests before creating section and allocating.
    CPPUNIT_ASSERT_EQUAL(PylithInt(0), field.chartSize());
    CPPUNIT_ASSERT_EQUAL(PylithInt(0), field.sectionSize());

    field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
    field.allocate();
    CPPUNIT_ASSERT_EQUAL(_data->numVertices, field.chartSize());
    CPPUNIT_ASSERT_EQUAL(_data->numVertices*fiberDim, field.sectionSize());

    PYLITH_METHOD_END;
} // testSectionAccessors

// ----------------------------------------------------------------------
// Test localVector(), globalVector().
void
pylith::topology::TestFieldMesh::testVectorAccessors(void)
{ // testVectorAccessors
    PYLITH_METHOD_BEGIN;

    Mesh mesh;
    _initializeMesh(&mesh);

    const int fiberDim = 2;
    const std::string& label = "field A";
    Field field(mesh);
    field.label(label.c_str());

    // Tests before creating section and allocating.
    CPPUNIT_ASSERT(!field.localVector());
    CPPUNIT_ASSERT(!field.globalVector());

    field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
    field.allocate();
    CPPUNIT_ASSERT(field.localVector());
    CPPUNIT_ASSERT(field.globalVector());

    PYLITH_METHOD_END;
} // testVectorAccessors

// ----------------------------------------------------------------------
// Test newSection(points), newSection(domain), newSetion(field).
void
pylith::topology::TestFieldMesh::testNewSection(void)
{ // testNewSection
    PYLITH_METHOD_BEGIN;

    Mesh mesh;
    _initializeMesh(&mesh);

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    const int fiberDim = 2;
    const std::string& label = "field A";

    { // newSection(points)
        const int npts = (vEnd-vStart) / 2;
        int_array pointsIn(npts);
        int_array pointsOut(vEnd-vStart - npts);
        int count = 0;
        size_t iIn = 0;
        size_t iOut = 0;
        for (PetscInt v = vStart; v < vEnd; ++v) {
            if (count % 2  == 0)
                pointsIn[iIn++] = v;
            else
                pointsOut[iOut++] = v;
            ++count;
        } // for
        CPPUNIT_ASSERT_EQUAL(iIn, pointsIn.size());
        CPPUNIT_ASSERT_EQUAL(iOut, pointsOut.size());

        Field field(mesh);
        field.label(label.c_str());
        field.newSection(pointsIn, fiberDim);
        field.allocate();

        // Points in array should have a fiber dimension of fiberDim.
        VecVisitorMesh fieldVisitor(field);
        for (size_t i = 0; i < pointsIn.size(); ++i) {
            CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(pointsIn[i]));
        } // for

        // Points not int array should have a fiber dimension of zero.
        PetscInt pStart, pEnd;
        PetscSection section = field.localSection(); CPPUNIT_ASSERT(section);
        PetscErrorCode err = PetscSectionGetChart(section, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        for (size_t i = 0; i < pointsOut.size(); ++i) {
            if (pointsOut[i] >= pStart && pointsOut[i] < pEnd) {
                CPPUNIT_ASSERT_EQUAL(0, fieldVisitor.sectionDof(pointsOut[i]));
            } // if
        } // for
        const char *name = NULL;
        PetscVec vec = field.localVector(); CPPUNIT_ASSERT(vec);
        err = PetscObjectGetName((PetscObject) vec, &name); PYLITH_CHECK_ERROR(err);
        CPPUNIT_ASSERT_EQUAL(label, std::string(name));
    } // newSection(points)

    Field field(mesh);
    { // newSection(domain)
        const std::string& label = "field A";
        field.label(label.c_str());
        field.newSection(Field::VERTICES_FIELD, fiberDim);
        field.allocate();

        PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
        Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
        const PetscInt vStart = depthStratum.begin();
        const PetscInt vEnd = depthStratum.end();

        VecVisitorMesh fieldVisitor(field);
        for(PetscInt v = vStart; v < vEnd; ++v) {
            CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
        } // for

        const char *name = NULL;
        PetscVec vec = field.localVector(); CPPUNIT_ASSERT(vec);
        PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name); PYLITH_CHECK_ERROR(err);
        CPPUNIT_ASSERT_EQUAL(label, std::string(name));
    } // newSection(domain)

    { // newSection(Field)
        const int fiberDim2 = 5;
        const std::string& label = "field B";
        Field field2(mesh);
        field2.label(label.c_str());
        field2.newSection(field, fiberDim2);
        field2.allocate();

        VecVisitorMesh field2Visitor(field2);
        for (PetscInt v = vStart; v < vEnd; ++v) {
            CPPUNIT_ASSERT_EQUAL(fiberDim2, field2Visitor.sectionDof(v));
        } // for

        const char *name = NULL;
        PetscVec vec = field2.localVector(); CPPUNIT_ASSERT(vec);
        PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name); PYLITH_CHECK_ERROR(err);
        CPPUNIT_ASSERT_EQUAL(label, std::string(name));
    } // newSection(Field)


    PYLITH_METHOD_END;
} // testNewSection

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
    _initializeMesh(&mesh);

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    PetscErrorCode err = 0;

    // Create field with atlas to use to create new field
    Field fieldSrc(mesh);
    { // Setup source field
        fieldSrc.newSection(Field::VERTICES_FIELD, fiberDim);
        PetscSection section = fieldSrc.localSection(); CPPUNIT_ASSERT(section);
        for(PetscInt v = vStart, iV = 0; v < vEnd; ++v) {
            err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]); PYLITH_CHECK_ERROR(err);
        } // for
        fieldSrc.allocate();

        int index = 0;
        for(PetscInt v = vStart, iV = 0; v < vEnd; ++v, index += nconstraints[iV++]) {
            err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[index]); PYLITH_CHECK_ERROR(
                err);
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
    CPPUNIT_ASSERT(section); CPPUNIT_ASSERT(vec);
    for(PetscInt v = vStart, iV = 0; v < vEnd; ++v) {
        PetscInt dof, cdof;
        err = PetscSectionGetDof(section, v, &dof); PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetConstraintDof(section, v, &cdof); PYLITH_CHECK_ERROR(err);
        CPPUNIT_ASSERT_EQUAL(fiberDim, dof);
        CPPUNIT_ASSERT_EQUAL(nconstraints[iV++], cdof);
    } // for

    // Verify vector scatters were also copied.
    CPPUNIT_ASSERT_EQUAL(fieldSrc._scatters[""].dm,  field._scatters[""].dm);
    CPPUNIT_ASSERT_EQUAL(fieldSrc._scatters["A"].dm, field._scatters["A"].dm);
    const char *name = NULL;
    err = PetscObjectGetName((PetscObject) vec, &name); PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(label, std::string(name));

    PYLITH_METHOD_END;
} // testCloneSection

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::testSubfields(void)
{ // testSubfields
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad Test not implemented.", false);

    PYLITH_METHOD_END;
} /// testSubfields

// ----------------------------------------------------------------------
// Test clear().
void
pylith::topology::TestFieldMesh::testClear(void)
{ // testClear
    PYLITH_METHOD_BEGIN;

    Mesh mesh(_data->cellDim);
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
    _initializeMesh(&mesh);

    Field field(mesh);
    field.scale(scale);
    field.newSection(Field::VERTICES_FIELD, fiberDim);
    field.allocate();

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
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
// Test zero(). zeroAll().
void
pylith::topology::TestFieldMesh::testZero(void)
{ // testZero
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
    _initializeMesh(&mesh);

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    // Create field and set constraint sizes
    Field field(mesh);
    field.scale(scale);
    field.newSection(Field::VERTICES_FIELD, fiberDim);
    PetscSection section = field.localSection(); CPPUNIT_ASSERT(section);
    PetscErrorCode err = 0;
    PetscInt index = 0;
    for(PetscInt v = vStart, iV=0; v < vEnd; ++v) {
        err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]); PYLITH_CHECK_ERROR(err);
    } // for
    field.allocate();
    PetscVec vec = field.localVector(); CPPUNIT_ASSERT(vec);
    for(PetscInt v = vStart, iV = 0; v < vEnd; ++v, index += nconstraints[iV++]) {
        err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[index]); PYLITH_CHECK_ERROR(err);
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
} // testZero

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
    _initializeMesh(&mesh);

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    Field field(mesh);
    field.scale(scale);
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
    _initializeMesh(&mesh);

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
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
    field.scale(scale);
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
    const char* fieldA = "temperature";
    const int fiberDimB = 2;
    const char* fieldB = "displacement";
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
    _initializeMesh(&mesh);

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    Field fieldSrc(mesh);
    { // Setup source field
        const char* temperatureComponents[1] = {"temperature"};
        const char* displacementComponents[2] = {"displacement_x", "displacement_y"};
        fieldSrc.label("solution");
        fieldSrc.subfieldAdd(fieldA, temperatureComponents, 1, Field::SCALAR, 1, 1, scale, true);
        fieldSrc.subfieldAdd(fieldB, displacementComponents, 2, Field::VECTOR, 1, 1, scale, true);
        fieldSrc.subfieldsSetup();
#if 0
        fieldSrc.newSection(Field::VERTICES_FIELD, fiberDim);
        fieldSrc.subfieldSetDof("one", Field::VERTICES_FIELD, fiberDimA);
        fieldSrc.subfieldSetDof("two", Field::VERTICES_FIELD, fiberDimB);
#endif
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
        field.copySubfield(fieldSrc, fieldB);

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
        field.copySubfield(fieldSrc, fieldA);

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
    _initializeMesh(&mesh);

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    Field fieldSrc(mesh);
    fieldSrc.scale(scale);
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
    field.scale(scale);
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
    _initializeMesh(&mesh);

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
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
    _initializeMesh(&mesh);

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    Field field(mesh);
    field.scale(scale);
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
pylith::topology::TestFieldMesh::testScatter(void)
{ // testScatter
    PYLITH_METHOD_BEGIN;

    const char* context = "abc";
    const int fiberDim = 3;
    const PylithScalar valuesE[4*3] = {
        1.1, 2.2, 3.3,
        1.2, 2.3, 3.4,
        1.3, 2.4, 3.5,
        1.4, 2.5, 3.6,
    };

    Mesh mesh;
    _initializeMesh(&mesh);

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    Field field(mesh);
    const std::string& label = "field A";
    field.label(label.c_str());
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

    { // Test createScatter(), scatterVector().
        CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());

        field.createScatter(mesh);
        CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
        const Field::ScatterInfo& sinfo = field._getScatter("");
        CPPUNIT_ASSERT(sinfo.dm);
        const PetscVec scatterVector = field.scatterVector();
        CPPUNIT_ASSERT_EQUAL(sinfo.vector, scatterVector);

        const int sizeE = (vEnd-vStart) * fiberDim;
        PylithInt size = 0;
        VecGetSize(scatterVector, &size);
        CPPUNIT_ASSERT_EQUAL(sizeE, size);

        // Check vector name
        const char* vecname = 0;
        PetscObjectGetName((PetscObject)scatterVector, &vecname);
        CPPUNIT_ASSERT_EQUAL(label, std::string(vecname));

        // Make sure we can do multiple calls to createScatter().
        field.createScatter(mesh);
        CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
    } // Test createScatter(), scatterVec().

    { // Test createScatterWithBC(), scatterLocalToContext(), scatterVectorToLocal().
        field.createScatterWithBC(mesh, context);
        field.scatterLocalToContext(context);

        const PetscVec vec = field.scatterVector(context); CPPUNIT_ASSERT(vec);
        PylithInt size = 0;
        PetscErrorCode err = VecGetSize(vec, &size); PYLITH_CHECK_ERROR(err);
        PylithScalar* valuesVec = NULL;
        err = VecGetArray(vec, &valuesVec); PYLITH_CHECK_ERROR(err);

        const PylithScalar tolerance = 1.0e-06;
        PylithInt sizeE = (vEnd-vStart) * fiberDim;
        CPPUNIT_ASSERT_EQUAL(sizeE, size);
        for (PylithInt i=0; i < sizeE; ++i) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i], valuesVec[i], tolerance);
        } // for
        err = VecRestoreArray(vec, &valuesVec); PYLITH_CHECK_ERROR(err);

        const PylithScalar scale = 0.25;
        err = VecScale(vec, scale);
        field.scatterVectorToLocal(vec, context);

        const PetscVec localVec = field.localVector(); CPPUNIT_ASSERT(localVec);
        size = 0;
        err = VecGetSize(localVec, &size); PYLITH_CHECK_ERROR(err);
        valuesVec = NULL;
        err = VecGetArray(localVec, &valuesVec); PYLITH_CHECK_ERROR(err);

        sizeE = (vEnd-vStart) * fiberDim;
        CPPUNIT_ASSERT_EQUAL(sizeE, size);
        for (PylithInt i=0; i < sizeE; ++i) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i]*scale, valuesVec[i], tolerance);
        } // for
        err = VecRestoreArray(localVec, &valuesVec); PYLITH_CHECK_ERROR(err);
    } // Test createScatterWithBC(), scatterLocalToContext(), scatterVectorToLocal().


    // Make sure scatters are replicated with cloneSection().
    Field field2(mesh);
    field2.cloneSection(field);
    CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

    const Field::ScatterInfo& sinfo2 = field2._getScatter("");
    CPPUNIT_ASSERT(sinfo2.dm);
    CPPUNIT_ASSERT(sinfo2.vector);

    const Field::ScatterInfo& sinfo2B = field2._getScatter(context);
    CPPUNIT_ASSERT(sinfo2B.dm);
    CPPUNIT_ASSERT(sinfo2B.vector);


    PYLITH_METHOD_END;
} // testScatter

#if 0 // :TODO: @brad Remove, obsolete.
// ----------------------------------------------------------------------
// Test splitDefault().
void
pylith::topology::TestFieldMesh::testSplitDefault(void)
{ // testSplitDefault
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    const PetscInt spaceDim =_data->cellDim;
    const PetscInt numFields = spaceDim;
    const PetscInt nconstraints[4] = { 1, 2, 0, 1 };
    const PetscInt constraints[4] = {
        1,     // 0
        0, 1,  // 1
               // 2
        0,     // 3
    };

    Mesh mesh;
    _initializeMesh(&mesh);

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
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
        PetscSection section = fieldSrc.localSection(); CPPUNIT_ASSERT(section);
        PetscInt iV = 0, iC = 0;
        for(PetscInt v = vStart; v < vEnd; ++v) {
            const int nconstraintsVertex = nconstraints[iV];
            err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]); PYLITH_CHECK_ERROR(err);
            for(PetscInt iConstraint = 0; iConstraint < nconstraintsVertex; ++iConstraint) {
                const int field = constraints[iC++];
                err = PetscSectionAddFieldConstraintDof(section, v, field, 1); PYLITH_CHECK_ERROR(err);
            } // for
        } // for
        fieldSrc.allocate();

        iV = 0; iC = 0;
        PetscInt zero = 0;
        for(PetscInt v = vStart; v < vEnd; ++v) {
            const int nconstraintsVertex = nconstraints[iV++];
            if (nconstraintsVertex > 0) {
                err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[iC]); PYLITH_CHECK_ERROR(
                    err);
            } // if
            for(PetscInt iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
                const PetscInt field = constraints[iC++];
                err = PetscSectionSetFieldConstraintIndices(section, v, field, &zero); PYLITH_CHECK_ERROR(err);
            } // for
        } // for
    } // Setup source field

    PetscSection section = fieldSrc.localSection(); CPPUNIT_ASSERT(section);
    PetscInt numSectionFields;
    err = PetscSectionGetNumFields(section, &numSectionFields); PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(numFields, numSectionFields);

    for(PetscInt f = 0; f < numFields; ++f) {
        PetscInt iV = 0, iC = 0;

        for(PetscInt v = vStart; v < vEnd; ++v, ++iV) {
            PetscInt fdof, fcdof;
            PetscBool isConstrained = PETSC_FALSE;
            const PetscInt nconstraintsVertex = nconstraints[iV];

            err = PetscSectionGetFieldDof(section, v, f, &fdof); PYLITH_CHECK_ERROR(err);
            err = PetscSectionGetFieldConstraintDof(section, v, f, &fcdof); PYLITH_CHECK_ERROR(err);
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
    CPPUNIT_ASSERT(_data);

    const PetscInt spaceDim = _data->cellDim;
    const PetscInt numFields = spaceDim;
    const PetscInt nconstraints[4] = { 1, 2, 0, 1 };
    const PetscInt constraints[4] = {
        1,     // 0
        0, 1,  // 1
               // 2
        0,     // 3
    };

    Mesh mesh;
    _initializeMesh(&mesh);

    PetscDM dmMesh = mesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
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
            err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]); PYLITH_CHECK_ERROR(err);
            for(PetscInt iConstraint = 0; iConstraint < nconstraintsVertex; ++iConstraint) {
                const int field = constraints[iC++];
                err = PetscSectionAddFieldConstraintDof(section, v, field, 1); PYLITH_CHECK_ERROR(err);
            } // for
        } // for
        fieldSrc.allocate();

        iV = 0; iC = 0;
        PetscInt zero = 0;
        for(PetscInt v = vStart; v < vEnd; ++v) {
            const int nconstraintsVertex = nconstraints[iV++];
            if (nconstraintsVertex > 0) {
                err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[iC]); PYLITH_CHECK_ERROR(
                    err);
            } // if
            for(PetscInt iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
                const PetscInt field = constraints[iC++];
                err = PetscSectionSetFieldConstraintIndices(section, v, field, &zero); PYLITH_CHECK_ERROR(err);
            } // for
        } // for
    } // Setup source field

    Field field(mesh);
    field.cloneSection(fieldSrc);

    PetscSection section = field.localSection(); CPPUNIT_ASSERT(section);
    PetscInt numSectionFields;
    err = PetscSectionGetNumFields(section, &numSectionFields); PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(numFields, numSectionFields);

    for(PetscInt f = 0; f < numFields; ++f) {
        PetscInt iV = 0, iC = 0;

        for(PetscInt v = vStart; v < vEnd; ++v, ++iV) {
            PetscInt fdof, fcdof;
            PetscBool isConstrained = PETSC_FALSE;
            const PetscInt nconstraintsVertex = nconstraints[iV];

            err = PetscSectionGetFieldDof(section, v, f, &fdof); PYLITH_CHECK_ERROR(err);
            err = PetscSectionGetFieldConstraintDof(section, v, f, &fcdof); PYLITH_CHECK_ERROR(err);
            CPPUNIT_ASSERT_EQUAL(1, fdof);
            for(PetscInt iConstraint = 0; iConstraint < nconstraintsVertex; ++iConstraint)
                if (constraints[iC++] == f) isConstrained = PETSC_TRUE;
            const PetscInt constraintDimE = (!isConstrained) ? 0 : 1;
            CPPUNIT_ASSERT_EQUAL(constraintDimE, fcdof);
        } // for
    } // for

    PYLITH_METHOD_END;
} // testCloneSectionSplit
#endif

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_initializeMesh(Mesh* mesh)
{ // _initializeMesh
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(mesh);

    const int cellDim = _data->cellDim;
    const int numCells = _data->numCells;
    const int numVertices = _data->numVertices;
    const int numCorners = _data->numCorners;
    const int spaceDim = _data->cellDim;

    PetscInt size = numVertices * spaceDim;
    scalar_array coordinates(size);
    for (PetscInt i=0; i < size; ++i) {
        coordinates[i] = _data->coordinates[i];
    }   // for

    size = numCells * numCorners;
    int_array cells(size);
    for (PetscInt i=0; i < size; ++i) {
        cells[i] = _data->cells[i];
    }   // for

    pylith::meshio::MeshBuilder::buildMesh(mesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners,
                                           cellDim);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(spaceDim);
    cs.initialize();
    mesh->coordsys(&cs);

    PYLITH_METHOD_END;
} // _initializeMesh


// ----------------------------------------------------------------------
// Constructor
pylith::topology::TestFieldMesh_Data::TestFieldMesh_Data(void) :
    cellDim(0),
    numVertices(0),
    numCells(0),
    numCorners(0),
    cells(NULL),
    coordinates(NULL)
{   // constructor
}   // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::topology::TestFieldMesh_Data::~TestFieldMesh_Data(void)
{   // destructor
}   // destructor


// End of file
