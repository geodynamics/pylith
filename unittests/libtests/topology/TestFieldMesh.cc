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
    CPPUNIT_ASSERT_EQUAL(_data->numVertices * fiberDim, size); // :TODO: @brad Fix this. (account for constraints)

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
    } // newSection(points)

    Field field(*_mesh);
    { // newSection(domain)
        const std::string& label = "field A";
        field.label(label.c_str());
        field.newSection(Field::VERTICES_FIELD, fiberDim);
        field.allocate();

        PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
        Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
        const PetscInt vStart = depthStratum.begin();
        const PetscInt vEnd = depthStratum.end();

        VecVisitorMesh fieldVisitor(field);
        for(PetscInt v = vStart; v < vEnd; ++v) {
            CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
        } // for

        const char *name = NULL;
        PetscVec vec = field.localVector(); CPPUNIT_ASSERT(vec);
        PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(label, std::string(name));
    } // newSection(domain)

    { // newSection(Field)
        const int fiberDim2 = 5;
        const std::string& label = "field B";
        Field field2(*_mesh);
        field2.label(label.c_str());
        field2.newSection(field, fiberDim2);
        field2.allocate();

        VecVisitorMesh field2Visitor(field2);
        for (PetscInt v = vStart; v < vEnd; ++v) {
            CPPUNIT_ASSERT_EQUAL(fiberDim2, field2Visitor.sectionDof(v));
        } // for

        const char *name = NULL;
        PetscVec vec = field2.localVector(); CPPUNIT_ASSERT(vec);
        PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name); CPPUNIT_ASSERT(!err);
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

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);

    PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

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
        PetscInt dof, cdof;
        err = PetscSectionGetDof(section, v, &dof); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(fiberDim, dof);
        err = PetscSectionGetConstraintDof(section, v, &cdof); CPPUNIT_ASSERT(!err);
        const PylithInt numConstraints = _data->subfieldANumConstraints[iV] + _data->subfieldBNumConstraints[iV];
        CPPUNIT_ASSERT_EQUAL(numConstraints, cdof);
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
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    { // Test with preallocated field
        Field field(*_mesh);
        field.newSection(Field::VERTICES_FIELD, _data->subfieldBNumComponents);
        field.allocate();
        field.copySubfield(*_field, _data->subfieldBName);

        VecVisitorMesh fieldVisitor(field);
        const PetscScalar* fieldArray = fieldVisitor.localArray();
        const PylithScalar tolerance = 1.0e-6;
        const PylithInt fiberDimE = _data->subfieldBNumComponents;
        for (PetscInt v = vStart, i = 0; v < vEnd; ++v) {
            const PetscInt off = fieldVisitor.sectionOffset(v);
            CPPUNIT_ASSERT_EQUAL(fiberDimE, fieldVisitor.sectionDof(v));
            for (PetscInt d = 0; d < fiberDimE; ++d) {
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
        for (PetscInt v = vStart, i = 0; v < vEnd; ++v) {
            const PetscInt off = fieldVisitor.sectionOffset(v);
            CPPUNIT_ASSERT_EQUAL(fiberDimE, fieldVisitor.sectionDof(v));
            for (PetscInt d = 0; d < fiberDimE; ++d) {
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
    _field->dimensionalizeOkay(true);
    _field->dimensionalize();

    PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    // Expect no change for this serial test.
    VecVisitorMesh fieldVisitor(*_field);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    const PylithInt fiberDim = _data->subfieldANumComponents + _data->subfieldBNumComponents;
    const PylithReal tolerance = 1.0e-6;
    for (PetscInt v=vStart, indexA=0, indexB=0; v < vEnd; ++v) {
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));

        // Check field A values.
        const PetscInt offA = fieldVisitor.sectionOffset(v);
        const PylithReal scaleA = _data->subfieldAScale;
        for (PetscInt d = 0; d < _data->subfieldANumComponents; ++d) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->subfieldAValues[indexA++]*scaleA, fieldArray[offA+d], tolerance);
        } // for

        // Check field B values.
        const PetscInt offB = offA + _data->subfieldANumComponents;
        const PylithReal scaleB = _data->subfieldBScale;
        for (PetscInt d = 0; d < _data->subfieldBNumComponents; ++d) {
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
        _field->scatterVectorToLocal(vec, context);
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

    delete _mesh; _mesh = new Mesh; CPPUNIT_ASSERT(_mesh);
    pylith::meshio::MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners,
                                           cellDim);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(spaceDim);
    cs.initialize();
    _mesh->coordsys(&cs);

    delete _field; _field = new Field(*_mesh);
    _field->label("solution");
    _field->subfieldAdd(_data->subfieldAName, _data->subfieldAComponents, _data->subfieldANumComponents,
                        _data->subfieldAType, _data->subfieldABasisOrder, _data->subfieldAQuadOrder,
                        _data->subfieldAScale, true);
    _field->subfieldAdd(_data->subfieldBName, _data->subfieldBComponents, _data->subfieldBNumComponents,
                        _data->subfieldBType, _data->subfieldBBasisOrder,_data->subfieldBQuadOrder,
                        _data->subfieldBScale, true);
    _field->subfieldsSetup();

    PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();
    PetscErrorCode err;

    // Set number of constraints
#if 0
    PetscSection section = _field->localSection(); CPPUNIT_ASSERT(section);
    for (PylithInt v=vStart, iVertex=0; v < vEnd; ++v) {
        const PylithInt numConstraintsA = _data->subfieldANumConstraints[iVertex];
        const PylithInt numConstraintsB = _data->subfieldBNumConstraints[iVertex];
        const PylithInt numConstraintsVertex = numConstraintsA + numConstraintsB;
        err = PetscSectionSetConstraintDof(section, v, numConstraintsVertex); CPPUNIT_ASSERT(!err);
        err = PetscSectionSetFieldConstraintDof(section, v, 0, numConstraintsA);
        err = PetscSectionSetFieldConstraintDof(section, v, 1, numConstraintsB);
    } // for
#endif

    // Allocate field.
    _field->allocate();
#if 0
    // Set constraint DOF.
    for (PylithInt v=vStart, iVertex=0, indexA=0, indexB=0; v < vEnd; ++v) {
        const PylithInt numConstraintsA = _data->subfieldANumConstraints[iVertex];
        const PylithInt numConstraintsB = _data->subfieldBNumConstraints[iVertex];
        const PylithInt numConstraintsVertex = numConstraintsA + numConstraintsB;
        int_array constraints(numConstraintsVertex);
        PylithInt iC = 0;
        for (PylithInt i=0; i < numConstraintsA; ++i) {
            constraints[iC++] = _data->subfieldAConstraints[indexA++];
        } // for
        for (PylithInt i=0; i < numConstraintsB; ++i) {
            constraints[iC++] = _data->subfieldBConstraints[indexB++];
        } // for
        err = PetscSectionSetConstraintIndices(section, v, &constraints[0]); CPPUNIT_ASSERT(!err);
    } // for
#endif

    // Populate with values.
    VecVisitorMesh fieldVisitor(*_field);
    const PylithInt fiberDim = _data->subfieldANumComponents + _data->subfieldBNumComponents;
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for (PylithInt v=vStart, indexA=0, indexB=0; v < vEnd; ++v) {
        // Set values for field A
        const PylithInt offA = fieldVisitor.sectionOffset(v);
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
        for (PetscInt d = 0; d < _data->subfieldANumComponents; ++d) {
            fieldArray[offA+d] = _data->subfieldAValues[indexA++];
        } // for
          // Set values for field B
        const PylithInt offB = offA + _data->subfieldANumComponents;
        for (PetscInt d = 0; d < _data->subfieldBNumComponents; ++d) {
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
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    const PylithInt fiberDim = _data->subfieldANumComponents + _data->subfieldBNumComponents;
    const PylithReal tolerance = 1.0e-6;
    for (PetscInt v=vStart, index=0; v < vEnd; ++v) {
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
    coordinates(NULL)
{   // constructor
}   // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::topology::TestFieldMesh_Data::~TestFieldMesh_Data(void)
{   // destructor
}   // destructor


// End of file
