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
// Copyright (c) 2010-2017 University of California, Davis
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
    _mesh = NULL;
    _field = NULL;

    PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::topology::TestFieldMesh::tearDown(void)
{ // tearDown
    PYLITH_METHOD_BEGIN;

    delete _data; _data = NULL;
    delete _mesh; _mesh = NULL;
    delete _field; _field = NULL;

    PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldMesh::testConstructor(void)
{ // testConstructor
    PYLITH_METHOD_BEGIN;

    Mesh meshA;
    Field field(meshA);

    PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test mesh().
void
pylith::topology::TestFieldMesh::testMesh(void)
{ // testMesh
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    _initialize();
    CPPUNIT_ASSERT(_field);

    CPPUNIT_ASSERT_EQUAL(_data->cellDim, _field->mesh().dimension());

    PYLITH_METHOD_END;
} // testMesh

// ----------------------------------------------------------------------
// Test label(), vectorFieldType(), scale(), addDimensionOkay(), spaceDim().
void
pylith::topology::TestFieldMesh::testGeneralAccessors(void)
{ // testGeneralAccessors
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    _initialize();
    CPPUNIT_ASSERT(_field);

    // Test label()
    const std::string label = "velocity";
    _field->label(label.c_str());
    CPPUNIT_ASSERT_EQUAL(label, std::string(_field->label()));

    // Test vectorFieldType()
    _field->vectorFieldType(FieldBase::SCALAR);
    CPPUNIT_ASSERT_EQUAL(FieldBase::SCALAR, _field->_metadata.vectorFieldType);

    // Test scale()
    const PylithScalar scale = 2.0;
    _field->scale(scale);
    const PylithScalar tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(scale, _field->_metadata.scale, tolerance);

    // Test addDimensionOkay()
    CPPUNIT_ASSERT_EQUAL(false, _field->_metadata.dimsOkay);
    _field->dimensionalizeOkay(true);
    CPPUNIT_ASSERT_EQUAL(true, _field->_metadata.dimsOkay);

    // Test spaceDim()
    CPPUNIT_ASSERT_EQUAL(_data->cellDim, _field->spaceDim());

    PYLITH_METHOD_END;
} // testGeneralAccessors

// ----------------------------------------------------------------------
// Test chartSize(), sectionSize(), localSection(), globalSection().
void
pylith::topology::TestFieldMesh::testSectionAccessors(void)
{ // testSectionAccessors
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    _initialize();
    CPPUNIT_ASSERT(_field);

    CPPUNIT_ASSERT(_field->chartSize() > 0); // vertices + edges + faces + cells
    const PylithInt fiberDim = _data->subfieldANumComponents + _data->subfieldBNumComponents;
    CPPUNIT_ASSERT_EQUAL(_data->numVertices*fiberDim, _field->sectionSize());

    PYLITH_METHOD_END;
} // testSectionAccessors

// ----------------------------------------------------------------------
// Test localVector(), globalVector().
void
pylith::topology::TestFieldMesh::testVectorAccessors(void)
{ // testVectorAccessors
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_field);

    PetscErrorCode err;
    PylithInt size = 0;

    const PetscVec& localVec = _field->localVector();
    err = VecGetSize(localVec, &size); CPPUNIT_ASSERT(!err);
    const PylithInt fiberDim = _data->subfieldANumComponents + _data->subfieldBNumComponents;
    CPPUNIT_ASSERT_EQUAL(_data->numVertices * fiberDim, size);

    const PetscVec& globalVec = _field->globalVector();
    err = VecGetSize(globalVec, &size); CPPUNIT_ASSERT(!err);
    // Count total number of expected constraints.
    PylithInt numConstraints = _data->bcANumVertices * _data->bcANumConstrainedDOF +
                               _data->bcBNumVertices * _data->bcBNumConstrainedDOF;
    CPPUNIT_ASSERT_EQUAL(_data->numVertices * fiberDim - numConstraints, size); // :TODO: @brad Fix this. (account for constraints)

    PYLITH_METHOD_END;
} // testVectorAccessors

// ----------------------------------------------------------------------
// Test newSection(points), newSection(domain), newSetion(field).
void
pylith::topology::TestFieldMesh::testNewSection(void)
{ // testNewSection
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);

    PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PylithInt vStart = depthStratum.begin();
    const PylithInt vEnd = depthStratum.end();

    const int fiberDim = 2;
    const std::string& label = "field A";
    { // newSection(points)
        const int npts = (vEnd-vStart) / 2;
        int_array pointsIn(npts);
        int_array pointsOut(vEnd-vStart - npts);
        int count = 0;
        size_t iIn = 0;
        size_t iOut = 0;
        for (PylithInt v = vStart; v < vEnd; ++v) {
            if (count % 2  == 0) {
                pointsIn[iIn++] = v;
            } else {
                pointsOut[iOut++] = v;
            } // if/else
            ++count;
        } // for
        CPPUNIT_ASSERT_EQUAL(iIn, pointsIn.size());
        CPPUNIT_ASSERT_EQUAL(iOut, pointsOut.size());

        Field field(*_mesh);
        field.label(label.c_str());
        field.newSection(&pointsIn[0], pointsIn.size(), fiberDim);
        field.allocate();

        // Points in array should have a fiber dimension of fiberDim.
        VecVisitorMesh fieldVisitor(field);
        for (size_t i = 0; i < pointsIn.size(); ++i) {
            CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(pointsIn[i]));
        } // for

        // Points not int array should have a fiber dimension of zero.
        PylithInt pStart, pEnd;
        PetscSection section = field.localSection(); CPPUNIT_ASSERT(section);
        PetscErrorCode err = PetscSectionGetChart(section, &pStart, &pEnd); CPPUNIT_ASSERT(!err);
        for (size_t i = 0; i < pointsOut.size(); ++i) {
            if (pointsOut[i] >= pStart && pointsOut[i] < pEnd) {
                CPPUNIT_ASSERT_EQUAL(0, fieldVisitor.sectionDof(pointsOut[i]));
            } // if
        } // for
        const char *name = NULL;
        PetscVec vec = field.localVector(); CPPUNIT_ASSERT(vec);
        err = PetscObjectGetName((PetscObject) vec, &name); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(label, std::string(name));

        // Test trapping with negative fiber dimension.
        Field fieldB(*_mesh);
        CPPUNIT_ASSERT_THROW(fieldB.newSection(pointsIn, -1), std::runtime_error);

        // Test creation with empty point list.
        pointsIn.resize(0);
        Field fieldC(*_mesh);
        fieldC.newSection(pointsIn, fiberDim);
        CPPUNIT_ASSERT_EQUAL(0, fieldC.chartSize());
        CPPUNIT_ASSERT_EQUAL(0, fieldC.sectionSize());
    } // newSection(points)

    Field field(*_mesh);
    { // newSection(VERTICES_FIELD)
        const std::string& label = "field A";
        field.label(label.c_str());
        field.newSection(Field::VERTICES_FIELD, fiberDim);
        field.allocate();

        PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
        Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
        const PylithInt vStart = depthStratum.begin();
        const PylithInt vEnd = depthStratum.end();

        VecVisitorMesh fieldVisitor(field);
        for(PylithInt v = vStart; v < vEnd; ++v) {
            CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
        } // for

        const char *name = NULL;
        PetscVec vec = field.localVector(); CPPUNIT_ASSERT(vec);
        PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(label, std::string(name));
    } // newSection(VERTICES_FIELD)

    { // newSection(CELLS_FIELD)
        Field fieldA(*_mesh);
        const std::string& label = "field B";
        fieldA.label(label.c_str());
        fieldA.newSection(Field::CELLS_FIELD, fiberDim);
        fieldA.allocate();

        PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
        Stratum heightStratum(dmMesh, Stratum::HEIGHT, 0);
        const PylithInt cStart = heightStratum.begin();
        const PylithInt cEnd = heightStratum.end();

        VecVisitorMesh fieldVisitor(fieldA);
        for(PylithInt c=cStart; c < cEnd; ++c) {
            CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(c));
        } // for

        const char *name = NULL;
        PetscVec vec = fieldA.localVector(); CPPUNIT_ASSERT(vec);
        PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(label, std::string(name));
    } // newSection(CELLS_FIELD)

    { // newSection(POINTS_FIELD)
        Field fieldA(*_mesh);
        const std::string& label = "field C";
        fieldA.label(label.c_str());
        fieldA.newSection(Field::POINTS_FIELD, fiberDim);
        fieldA.allocate();

        PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
        PylithInt pStart, pEnd;
        PetscErrorCode err = DMPlexGetChart(dmMesh, &pStart, &pEnd); CPPUNIT_ASSERT(!err);

        VecVisitorMesh fieldVisitor(fieldA);
        for(PylithInt p=pStart; p < pEnd; ++p) {
            CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(p));
        } // for

        const char *name = NULL;
        PetscVec vec = fieldA.localVector(); CPPUNIT_ASSERT(vec);
        err = PetscObjectGetName((PetscObject) vec, &name); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(label, std::string(name));
    } // newSection(POINTS_FIELD)

    { // newSection(Field)
        const int fiberDim2 = 5;
        const std::string& label = "field B";
        Field field2(*_mesh);
        field2.label(label.c_str());
        field2.newSection(field, fiberDim2);
        field2.allocate();

        VecVisitorMesh field2Visitor(field2);
        for (PylithInt v = vStart; v < vEnd; ++v) {
            CPPUNIT_ASSERT_EQUAL(fiberDim2, field2Visitor.sectionDof(v));
        } // for

        const char *name = NULL;
        PetscVec vec = field2.localVector(); CPPUNIT_ASSERT(vec);
        PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(label, std::string(name));

        // Test trapping with negative fiber dimension.
        Field field3(*_mesh);
        CPPUNIT_ASSERT_THROW(field3.newSection(field, -1), std::runtime_error);
    } // newSection(Field)


    PYLITH_METHOD_END;
} // testNewSection

// ----------------------------------------------------------------------
// Test cloneSection().
void
pylith::topology::TestFieldMesh::testCloneSection(void)
{ // testCloneSection
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);

    PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PylithInt vStart = depthStratum.begin();
    const PylithInt vEnd = depthStratum.end();

    PetscErrorCode err = 0;

    // Add two scatters, one with default context and one with given context.
    _field->createScatter(*_mesh);
    const char* context = "ABC";
    _field->createScatter(*_mesh, context);

    Field field(*_mesh);
    const std::string& label = "field A";
    field.label(label.c_str());
    field.cloneSection(*_field);
    PetscSection section = field.localSection();
    PetscVec vec = field.localVector();
    CPPUNIT_ASSERT(section); CPPUNIT_ASSERT(vec);
    const PylithInt fiberDim = _data->subfieldANumComponents + _data->subfieldBNumComponents;
    for (PylithInt v=vStart, iV=0; v < vEnd; ++v, ++iV) {
        PylithInt dof, cdof;
        err = PetscSectionGetDof(section, v, &dof); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(fiberDim, dof);

        // Count number of expected constraints on vertex.
        PylithInt numConstraintsE = 0;
        for (int i=0; i < _data->bcANumVertices; ++i) {
            const PylithInt vIndex = v - _data->numCells;
            if (_data->bcAVertices[i] == vIndex) {
                numConstraintsE += _data->bcANumConstrainedDOF;
                break;
            }
        } // for
        for (int i=0; i < _data->bcBNumVertices; ++i) {
            const PylithInt vIndex = v - _data->numCells;
            if (_data->bcBVertices[i] == vIndex) {
                numConstraintsE += _data->bcBNumConstrainedDOF;
                break;
            }
        } // for
        err = PetscSectionGetConstraintDof(section, v, &cdof); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(numConstraintsE, cdof);
    } // for

    // Verify vector scatters were also copied.
    CPPUNIT_ASSERT_EQUAL(_field->_scatters[""].dm,  field._scatters[""].dm);
    CPPUNIT_ASSERT_EQUAL(_field->_scatters[context].dm, field._scatters[context].dm);
    const char *name = NULL;
    err = PetscObjectGetName((PetscObject) vec, &name); CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(label, std::string(name));

    PYLITH_METHOD_END;
} // testCloneSection

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::testSubfieldAccessors(void)
{ // testSubfieldAccessors
    PYLITH_METHOD_BEGIN;
    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);

    // Subfields setup via subfieldAdd() and subfieldsSetup() in _initialize().

    // Test hasSubfield().
    CPPUNIT_ASSERT(_field->hasSubfield(_data->subfieldAName));
    CPPUNIT_ASSERT(_field->hasSubfield(_data->subfieldBName));
    CPPUNIT_ASSERT(!_field->hasSubfield("zyxwvut987654321"));

    // Test subfieldNames().
    const string_vector& names = _field->subfieldNames();
    CPPUNIT_ASSERT_EQUAL(size_t(2), names.size());
    CPPUNIT_ASSERT_EQUAL(std::string(_data->subfieldAName), names[0]);
    CPPUNIT_ASSERT_EQUAL(std::string(_data->subfieldBName), names[1]);

    { // Test subfieldInfo() for subfield A.
        const Field::SubfieldInfo& infoA = _field->subfieldInfo(_data->subfieldAName);
        CPPUNIT_ASSERT_EQUAL(0, infoA.index);
        CPPUNIT_ASSERT_EQUAL(_data->subfieldANumComponents, infoA.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(_data->subfieldAName), infoA.metadata.label);
        CPPUNIT_ASSERT_EQUAL(_data->subfieldAType, infoA.metadata.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(_data->subfieldAScale, infoA.metadata.scale);
        const string_vector& componentNames = infoA.metadata.componentNames;
        CPPUNIT_ASSERT_EQUAL(size_t(_data->subfieldANumComponents), componentNames.size());
        for (int i=0; i < _data->subfieldANumComponents; ++i) {
            CPPUNIT_ASSERT_EQUAL(std::string(_data->subfieldAComponents[i]), componentNames[i]);
        } // for
        CPPUNIT_ASSERT_EQUAL(_data->subfieldABasisOrder, infoA.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(_data->subfieldAQuadOrder, infoA.fe.quadOrder);
    } // Test subfieldInfo() for subfield A.

    { // Test subfieldInfo() for subfield B.
        const Field::SubfieldInfo& infoB = _field->subfieldInfo(_data->subfieldBName);
        CPPUNIT_ASSERT_EQUAL(_data->subfieldBNumComponents, infoB.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(_data->subfieldBName), infoB.metadata.label);
        CPPUNIT_ASSERT_EQUAL(_data->subfieldBType, infoB.metadata.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(_data->subfieldBScale, infoB.metadata.scale);
        const string_vector& componentNames = infoB.metadata.componentNames;
        CPPUNIT_ASSERT_EQUAL(size_t(_data->subfieldBNumComponents), componentNames.size());
        for (int i=0; i < _data->subfieldBNumComponents; ++i) {
            CPPUNIT_ASSERT_EQUAL(std::string(_data->subfieldBComponents[i]), componentNames[i]);
        } // for
        CPPUNIT_ASSERT_EQUAL(_data->subfieldBBasisOrder, infoB.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(_data->subfieldBQuadOrder, infoB.fe.quadOrder);
    } // Test subfieldInfo() for subfield B.

    CPPUNIT_ASSERT_THROW(_field->subfieldInfo("aabbccdd"), std::runtime_error);

    PYLITH_METHOD_END;
} /// testSubfieldAccessors

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

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);

    _checkValues(*_field);

    PYLITH_METHOD_END;
} // testAllocate

// ----------------------------------------------------------------------
// Test zeroLocal().
void
pylith::topology::TestFieldMesh::testZeroLocal(void)
{ // testZeroLocal
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);
    _field->zeroLocal();

    _checkValues(*_field, 0.0);

    PYLITH_METHOD_END;
} // testZeroLocal

// ----------------------------------------------------------------------
// Test complete().
void
pylith::topology::TestFieldMesh::testComplete(void)
{ // testComplete
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);
    _field->complete();

    // Expect no change for this serial test.
    _checkValues(*_field);

    PYLITH_METHOD_END;
} // testComplete

// ----------------------------------------------------------------------
// Test copy().
void
pylith::topology::TestFieldMesh::testCopy(void)
{ // testCopy
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);

    Field field(*_mesh);
    field.cloneSection(*_field);
    field.allocate();
    field.copy(*_field);

    // Expect no change for this serial test.
    _checkValues(field);

    // Test trapping when wrong size.
    Field field2(*_mesh);
    field2.allocate();
    CPPUNIT_ASSERT_THROW(field2.copy(*_field), std::runtime_error);

    PYLITH_METHOD_END;
} // testCopy

// ----------------------------------------------------------------------
// Test copySubfield().
void
pylith::topology::TestFieldMesh::testCopySubfield(void)
{ // testCopy
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);

    PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PylithInt vStart = depthStratum.begin();
    const PylithInt vEnd = depthStratum.end();

    { // Test with preallocated field
        Field field(*_mesh);
        field.newSection(Field::VERTICES_FIELD, _data->subfieldBNumComponents);
        field.allocate();
        field.copySubfield(*_field, _data->subfieldBName);

        VecVisitorMesh fieldVisitor(field);
        const PetscScalar* fieldArray = fieldVisitor.localArray();
        const PylithScalar tolerance = 1.0e-6;
        const PylithInt fiberDimE = _data->subfieldBNumComponents;
        for (PylithInt v = vStart, i = 0; v < vEnd; ++v) {
            const PylithInt off = fieldVisitor.sectionOffset(v);
            CPPUNIT_ASSERT_EQUAL(fiberDimE, fieldVisitor.sectionDof(v));
            for (PylithInt d = 0; d < fiberDimE; ++d) {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->subfieldBValues[i++], fieldArray[off+d], tolerance);
            } // for
        } // for
    } // Test with preallocated field

    { // Test with unallocated field
        Field field(*_mesh);
        field.copySubfield(*_field, _data->subfieldBName);

        VecVisitorMesh fieldVisitor(field);
        const PetscScalar* fieldArray = fieldVisitor.localArray();
        const PylithScalar tolerance = 1.0e-6;
        const PylithInt fiberDimE = _data->subfieldBNumComponents;
        for (PylithInt v = vStart, i = 0; v < vEnd; ++v) {
            const PylithInt off = fieldVisitor.sectionOffset(v);
            CPPUNIT_ASSERT_EQUAL(fiberDimE, fieldVisitor.sectionDof(v));
            for (PylithInt d = 0; d < fiberDimE; ++d) {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->subfieldBValues[i++], fieldArray[off+d], tolerance);
            } // for
        } // for
    } // Test with unallocated field

    PYLITH_METHOD_END;
} // testCopySubfield

// ----------------------------------------------------------------------
// Test dimensionalize().
void
pylith::topology::TestFieldMesh::testDimensionalize(void)
{ // testDimensionalize
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);

    // Default is not okay, so verify error trapping.
    CPPUNIT_ASSERT_THROW(_field->dimensionalize(), std::runtime_error);

    // Enable dimensionalization and then test result.
    _field->dimensionalizeOkay(true);
    _field->dimensionalize();

    PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PylithInt vStart = depthStratum.begin();
    const PylithInt vEnd = depthStratum.end();

    // Expect no change for this serial test.
    VecVisitorMesh fieldVisitor(*_field);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    const PylithInt fiberDim = _data->subfieldANumComponents + _data->subfieldBNumComponents;
    const PylithReal tolerance = 1.0e-6;
    for (PylithInt v=vStart, indexA=0, indexB=0; v < vEnd; ++v) {
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));

        // Check field A values.
        const PylithInt offA = fieldVisitor.sectionOffset(v);
        const PylithReal scaleA = _data->subfieldAScale;
        for (PylithInt d = 0; d < _data->subfieldANumComponents; ++d) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->subfieldAValues[indexA++]*scaleA, fieldArray[offA+d], tolerance);
        } // for

        // Check field B values.
        const PylithInt offB = offA + _data->subfieldANumComponents;
        const PylithReal scaleB = _data->subfieldBScale;
        for (PylithInt d = 0; d < _data->subfieldBNumComponents; ++d) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->subfieldBValues[indexB++]*scaleB, fieldArray[offB+d], tolerance);
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

    _initialize();
    CPPUNIT_ASSERT(_field);
    _field->view("Testing view");

    PYLITH_METHOD_END;
} // testView

// ----------------------------------------------------------------------
// Test createScatter().
void
pylith::topology::TestFieldMesh::testScatter(void)
{ // testScatter
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);

    { // Test createScatter(), scatterVector().
        CPPUNIT_ASSERT_EQUAL(size_t(0), _field->_scatters.size());

        _field->createScatter(*_mesh);
        CPPUNIT_ASSERT_EQUAL(size_t(1), _field->_scatters.size());
        const Field::ScatterInfo& sinfo = _field->_getScatter("");
        CPPUNIT_ASSERT(sinfo.dm);
        const PetscVec scatterVector = _field->scatterVector();
        CPPUNIT_ASSERT_EQUAL(sinfo.vector, scatterVector);
        // Check vector name
        const char* vecname = 0;
        PetscErrorCode err = PetscObjectGetName((PetscObject)scatterVector, &vecname); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(std::string("solution"), std::string(vecname));

        // Make sure we can do multiple calls to createScatter().
        _field->createScatter(*_mesh);
        CPPUNIT_ASSERT_EQUAL(size_t(1), _field->_scatters.size());
    } // Test createScatter(), scatterVec().

    const char* context = "ABC";
    { // Test createScatterWithBC(), scatterLocalToContext(), scatterVectorToLocal().
        _field->createScatterWithBC(*_mesh, context);
        _field->scatterLocalToContext(context);

        const PetscVec& vec = _field->scatterVector(context); CPPUNIT_ASSERT(vec);
        _checkValues(vec);

        const PylithScalar scale = 0.25;
        PetscErrorCode err = VecScale(vec, scale); CPPUNIT_ASSERT(!err);
        _field->scatterContextToLocal(context);
        const PetscVec localVec = _field->localVector(); CPPUNIT_ASSERT(localVec);
        _checkValues(localVec, scale);
    } // Test createScatterWithBC(), scatterLocalToContext(), scatterVectorToLocal().

    // Make sure scatters are replicated with cloneSection().
    Field field2(*_mesh);
    field2.cloneSection(*_field);
    CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

    const Field::ScatterInfo& sinfo2 = field2._getScatter("");
    CPPUNIT_ASSERT(sinfo2.dm);
    CPPUNIT_ASSERT(sinfo2.vector);

    const Field::ScatterInfo& sinfo2B = field2._getScatter(context);
    CPPUNIT_ASSERT(sinfo2B.dm);
    CPPUNIT_ASSERT(sinfo2B.vector);

    CPPUNIT_ASSERT_THROW(field2._getScatter("zzyyxx"), std::runtime_error);

    PYLITH_METHOD_END;
} // testScatter


// ----------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_initialize(void)
{ // _initialize
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    const int cellDim = _data->cellDim;
    const int numCells = _data->numCells;
    const int numVertices = _data->numVertices;
    const int numCorners = _data->numCorners;
    const int spaceDim = _data->cellDim;

    PylithInt size = numVertices * spaceDim;
    scalar_array coordinates(size);
    for (PylithInt i=0; i < size; ++i) {
        coordinates[i] = _data->coordinates[i];
    }   // for

    size = numCells * numCorners;
    int_array cells(size);
    for (PylithInt i=0; i < size; ++i) {
        cells[i] = _data->cells[i];
    }   // for

    delete _mesh; _mesh = new Mesh; CPPUNIT_ASSERT(_mesh);
    pylith::meshio::MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners,
                                           cellDim);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(spaceDim);
    cs.initialize();
    _mesh->coordsys(&cs);

    // Setup labels for constraints.
    PetscErrorCode err;
    for (PylithInt i=0; i < _data->bcANumVertices; ++i) {
        err = DMSetLabelValue(_mesh->dmMesh(), _data->bcALabel, numCells+_data->bcAVertices[i], _data->bcALabelId);
        CPPUNIT_ASSERT(!err);
    }   // for
    for (PylithInt i=0; i < _data->bcBNumVertices; ++i) {
        err = DMSetLabelValue(_mesh->dmMesh(), _data->bcBLabel, numCells+_data->bcBVertices[i], _data->bcBLabelId);
        CPPUNIT_ASSERT(!err);
    }   // for


    // Setup field
    delete _field; _field = new Field(*_mesh);
    _field->label("solution");
    _field->subfieldAdd(_data->subfieldAName, _data->subfieldAComponents, _data->subfieldANumComponents,
                        _data->subfieldAType, _data->subfieldABasisOrder, _data->subfieldAQuadOrder, true,
                        _data->subfieldAScale, NULL);
    _field->subfieldAdd(_data->subfieldBName, _data->subfieldBComponents, _data->subfieldBNumComponents,
                        _data->subfieldBType, _data->subfieldBBasisOrder, _data->subfieldBQuadOrder, true,
                        _data->subfieldBScale, NULL);
    _field->subfieldsSetup();

    err = DMAddBoundary(_field->dmMesh(), DM_BC_ESSENTIAL, "bcA", _data->bcALabel, 0, _data->bcANumConstrainedDOF,
                        _data->bcAConstrainedDOF, NULL, 1, &_data->bcALabelId, NULL); CPPUNIT_ASSERT(!err);
    err = DMAddBoundary(_field->dmMesh(), DM_BC_ESSENTIAL, "bcB", _data->bcBLabel, 0, _data->bcBNumConstrainedDOF,
                        _data->bcBConstrainedDOF, NULL, 1, &_data->bcBLabelId, NULL); CPPUNIT_ASSERT(!err);

    // Allocate field.
    _field->allocate();

    // Populate with values.
    PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PylithInt vStart = depthStratum.begin();
    const PylithInt vEnd = depthStratum.end();

    VecVisitorMesh fieldVisitor(*_field);
    const PylithInt fiberDim = _data->subfieldANumComponents + _data->subfieldBNumComponents;
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for (PylithInt v=vStart, indexA=0, indexB=0; v < vEnd; ++v) {
        // Set values for field A
        const PylithInt offA = fieldVisitor.sectionOffset(v);
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
        for (PylithInt d = 0; d < _data->subfieldANumComponents; ++d) {
            fieldArray[offA+d] = _data->subfieldAValues[indexA++];
        } // for
          // Set values for field B
        const PylithInt offB = offA + _data->subfieldANumComponents;
        for (PylithInt d = 0; d < _data->subfieldBNumComponents; ++d) {
            fieldArray[offB+d] = _data->subfieldBValues[indexB++];
        } // for
    } // for

    PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_checkValues(const Field& field,
                                              const PylithReal scale)
{ // _checkValues
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    // Create array of values in field from subfields.
    const int numVertices = _data->numVertices;
    const int fiberDimA = _data->subfieldANumComponents;
    const int fiberDimB = _data->subfieldBNumComponents;
    scalar_array valuesE(numVertices * (fiberDimA + fiberDimB));
    for (int iVertex=0, index=0, indexA=0, indexB=0; iVertex < numVertices; ++iVertex) {
        for (int d=0; d < fiberDimA; ++d) {
            valuesE[index++] = _data->subfieldAValues[indexA++];
        } // for
        for (int d=0; d < fiberDimB; ++d) {
            valuesE[index++] = _data->subfieldBValues[indexB++];
        } // for
    } // for

    PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PylithInt vStart = depthStratum.begin();
    const PylithInt vEnd = depthStratum.end();

    VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    const PylithInt fiberDim = _data->subfieldANumComponents + _data->subfieldBNumComponents;
    const PylithReal tolerance = 1.0e-6;
    for (PylithInt v=vStart, index=0; v < vEnd; ++v) {
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
        const PylithInt off = fieldVisitor.sectionOffset(v);

        for (PylithInt d=0; d < fiberDim; ++d) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[index++]*scale, fieldArray[off+d], tolerance);
        } // for
    } // for

    PYLITH_METHOD_END;
} // _checkValues

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_checkValues(const PetscVec& vec,
                                              const PylithReal scale)
{ // _checkValues
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    // Create array of values in field from subfields.
    const int numVertices = _data->numVertices;
    const int fiberDimA = _data->subfieldANumComponents;
    const int fiberDimB = _data->subfieldBNumComponents;
    scalar_array valuesE(numVertices * (fiberDimA + fiberDimB));
    for (int iVertex=0, index=0, indexA=0, indexB=0; iVertex < numVertices; ++iVertex) {
        for (int d=0; d < fiberDimA; ++d) {
            valuesE[index++] = _data->subfieldAValues[indexA++];
        } // for
        for (int d=0; d < fiberDimB; ++d) {
            valuesE[index++] = _data->subfieldBValues[indexB++];
        } // for
    } // for

    PetscErrorCode err;
    PylithInt size = 0;
    PylithScalar* vecArray = NULL;
    err = VecGetSize(vec, &size); CPPUNIT_ASSERT(!err);
    err = VecGetArray(vec, &vecArray); CPPUNIT_ASSERT(!err);

    const PylithInt sizeE = numVertices * (fiberDimA + fiberDimB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_EQUAL(sizeE, size);
    for (PylithInt i=0; i < sizeE; ++i) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i]*scale, vecArray[i], tolerance);
    } // for
    err = VecRestoreArray(vec, &vecArray); CPPUNIT_ASSERT(!err);

    PYLITH_METHOD_END;
} // _checkValues

// ----------------------------------------------------------------------
// Constructor
pylith::topology::TestFieldMesh_Data::TestFieldMesh_Data(void) :
    cellDim(0),
    numVertices(0),
    numCells(0),
    numCorners(0),
    cells(NULL),
    coordinates(NULL),

    subfieldAName(NULL),
    subfieldAType(Field::OTHER),
    subfieldAScale(0.0),
    subfieldANumComponents(0),
    subfieldAComponents(NULL),
    subfieldAValues(NULL),
    subfieldABasisOrder(0),
    subfieldAQuadOrder(0),
    bcALabel(NULL),
    bcALabelId(0),
    bcANumConstrainedDOF(0),
    bcAConstrainedDOF(NULL),
    bcANumVertices(0),
    bcAVertices(NULL),

    subfieldBName(NULL),
    subfieldBType(Field::OTHER),
    subfieldBScale(0.0),
    subfieldBNumComponents(0),
    subfieldBComponents(NULL),
    subfieldBValues(NULL),
    subfieldBBasisOrder(0),
    subfieldBQuadOrder(0),
    bcBLabel(NULL),
    bcBLabelId(0),
    bcBNumConstrainedDOF(0),
    bcBConstrainedDOF(NULL),
    bcBNumVertices(0),
    bcBVertices(NULL)
{   // constructor
}   // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::topology::TestFieldMesh_Data::~TestFieldMesh_Data(void)
{   // destructor
}   // destructor


// End of file
