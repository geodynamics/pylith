// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
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

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::topology::TestFieldMesh::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _data = new TestFieldMesh_Data;CPPUNIT_ASSERT(_data);
    _mesh = NULL;
    _field = NULL;

    PYLITH_METHOD_END;
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::topology::TestFieldMesh::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;
    delete _field;_field = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldMesh::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    Mesh meshA;
    Field field(meshA);

    PYLITH_METHOD_END;
} // testConstructor


// ---------------------------------------------------------------------------------------------------------------------
// Test cloneSection().
void
pylith::topology::TestFieldMesh::testCopyConstructor(void) {
    PYLITH_METHOD_BEGIN;
    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);

    PetscDM dmMesh = _mesh->getDM();CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PylithInt vStart = depthStratum.begin();
    const PylithInt vEnd = depthStratum.end();

    PetscErrorCode err = 0;

    const std::string& label = "field A";
    Field field(*_field);
    field.setLabel(label.c_str());

    const char *name = NULL;
    err = PetscObjectGetName((PetscObject)field.getDM(), &name);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(label, std::string(name));

    PetscSection section = field.getLocalSection();CPPUNIT_ASSERT(section);
    PetscVec vec = field.getLocalVector();CPPUNIT_ASSERT(vec);

    err = PetscObjectGetName((PetscObject) vec, &name);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(label, std::string(name));

    const PylithInt ndof = _data->descriptionA.numComponents + _data->descriptionB.numComponents;
    for (PylithInt v = vStart, iV = 0; v < vEnd; ++v, ++iV) {
        PylithInt dof, cdof;
        err = PetscSectionGetDof(section, v, &dof);CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(ndof, dof);

        // Count number of expected constraints on vertex.
        PylithInt numConstraintsE = 0;
        for (int i = 0; i < _data->bcANumVertices; ++i) {
            const PylithInt vIndex = v - _data->numCells;
            if (_data->bcAVertices[i] == vIndex) {
                numConstraintsE += _data->bcANumConstrainedDOF;
                break;
            }
        } // for
        for (int i = 0; i < _data->bcBNumVertices; ++i) {
            const PylithInt vIndex = v - _data->numCells;
            if (_data->bcBVertices[i] == vIndex) {
                numConstraintsE += _data->bcBNumConstrainedDOF;
                break;
            }
        } // for
        err = PetscSectionGetConstraintDof(section, v, &cdof);CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(numConstraintsE, cdof);
    } // for

    field.deallocate();

    PYLITH_METHOD_END;
} // testCopyConstructor


// ---------------------------------------------------------------------------------------------------------------------
// Test mesh().
void
pylith::topology::TestFieldMesh::testMesh(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    _initialize();
    CPPUNIT_ASSERT(_field);

    CPPUNIT_ASSERT_EQUAL(_data->cellDim, _field->getMesh().getDimension());

    PYLITH_METHOD_END;
} // testMesh


// ---------------------------------------------------------------------------------------------------------------------
// Test getLabel(), vectorFieldType(), scale(), addDimensionOkay(), getSpaceDim().
void
pylith::topology::TestFieldMesh::testGeneralAccessors(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    _initialize();
    CPPUNIT_ASSERT(_field);

    // Test getLabel()
    const std::string label = "velocity";
    _field->setLabel(label.c_str());
    CPPUNIT_ASSERT_EQUAL(label, std::string(_field->getLabel()));
    const char* name = NULL;
    PetscErrorCode err = 0;
    err = PetscObjectGetName((PetscObject)_field->getDM(), &name);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(label, std::string(name));

    // Test getSpaceDim()
    CPPUNIT_ASSERT_EQUAL(size_t(_data->cellDim), _field->getSpaceDim());

    PYLITH_METHOD_END;
} // testGeneralAccessors


// ---------------------------------------------------------------------------------------------------------------------
// Test chartSize(), getStorageSize(), getLocalSection(), globalSection().
void
pylith::topology::TestFieldMesh::testSectionAccessors(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    _initialize();
    CPPUNIT_ASSERT(_field);

    CPPUNIT_ASSERT(_field->getChartSize() > 0); // vertices + edges + faces + cells
    const PylithInt ndof = _data->descriptionA.numComponents + _data->descriptionB.numComponents;
    CPPUNIT_ASSERT_EQUAL(_data->numVertices*ndof, _field->getStorageSize());

    PYLITH_METHOD_END;
} // testSectionAccessors


// ---------------------------------------------------------------------------------------------------------------------
// Test getLocalVector(), getGlobalVector().
void
pylith::topology::TestFieldMesh::testVectorAccessors(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_field);

    PetscErrorCode err;
    const char* name = NULL;
    PylithInt size = 0;
    const PylithInt ndof =
        _data->numVertices * (_data->descriptionA.numComponents + _data->descriptionB.numComponents);
    const PylithInt ndofConstrained =
        _data->bcANumConstrainedDOF*_data->bcANumVertices + _data->bcBNumConstrainedDOF*_data->bcBNumVertices;

    const PetscVec& localVec = _field->getLocalVector();CPPUNIT_ASSERT(localVec);
    err = PetscObjectGetName((PetscObject)localVec, &name);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(std::string(_field->getLabel()), std::string(name));
    err = VecGetSize(localVec, &size);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(ndof, size);

    _field->createGlobalVector();
    const PetscVec& globalVec = _field->getGlobalVector();CPPUNIT_ASSERT(globalVec);
    _field->scatterLocalToVector(globalVec);
    err = PetscObjectGetName((PetscObject)globalVec, &name);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(std::string(_field->getLabel()), std::string(name));
    err = VecGetSize(globalVec, &size);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(ndof - ndofConstrained, size);

    _field->createOutputVector();
    _field->scatterLocalToOutput();
    const PetscVec& outputVec = _field->getOutputVector();CPPUNIT_ASSERT(outputVec);
    err = PetscObjectGetName((PetscObject)outputVec, &name);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(std::string(_field->getLabel()), std::string(name));
    err = VecGetSize(outputVec, &size);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(ndof, size);
    _checkValues(outputVec);

    PYLITH_METHOD_END;
} // testVectorAccessors


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::testSubfieldAccessors(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);

    // Subfields setup via subfieldAdd() and subfieldsSetup() in _initialize().

    // Test hasSubfield().
    CPPUNIT_ASSERT(_field->hasSubfield(_data->descriptionA.label.c_str()));
    CPPUNIT_ASSERT(_field->hasSubfield(_data->descriptionB.label.c_str()));
    CPPUNIT_ASSERT(!_field->hasSubfield("zyxwvut987654321"));

    // Test getSubfieldNames().
    const string_vector& names = _field->getSubfieldNames();
    CPPUNIT_ASSERT_EQUAL(size_t(2), names.size());
    CPPUNIT_ASSERT_EQUAL(_data->descriptionA.label, names[0]);
    CPPUNIT_ASSERT_EQUAL(_data->descriptionB.label, names[1]);

    { // Test getSubfieldInfo() for subfield A.
        const Field::SubfieldInfo& infoA = _field->getSubfieldInfo(_data->descriptionA.label.c_str());
        CPPUNIT_ASSERT_EQUAL(0, infoA.index);
        CPPUNIT_ASSERT_EQUAL(_data->descriptionA.numComponents, infoA.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(_data->descriptionA.label, infoA.description.label);
        CPPUNIT_ASSERT_EQUAL(_data->descriptionA.vectorFieldType, infoA.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(_data->descriptionA.scale, infoA.description.scale);
        const string_vector& componentNames = infoA.description.componentNames;
        CPPUNIT_ASSERT_EQUAL(_data->descriptionA.numComponents, componentNames.size());
        for (size_t i = 0; i < _data->descriptionA.numComponents; ++i) {
            CPPUNIT_ASSERT_EQUAL(_data->descriptionA.componentNames[i], componentNames[i]);
        } // for
        CPPUNIT_ASSERT_EQUAL(_data->discretizationA.basisOrder, infoA.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(_data->discretizationA.quadOrder, infoA.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(_data->discretizationA.feSpace, infoA.fe.feSpace);
        CPPUNIT_ASSERT_EQUAL(_data->discretizationA.isBasisContinuous, infoA.fe.isBasisContinuous);
    } // Test getSubfieldInfo() for subfield A.

    { // Test getSubfieldInfo() for subfield B.
        const Field::SubfieldInfo& infoB = _field->getSubfieldInfo(_data->descriptionB.label.c_str());
        CPPUNIT_ASSERT_EQUAL(1, infoB.index);
        CPPUNIT_ASSERT_EQUAL(_data->descriptionB.numComponents, infoB.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(_data->descriptionB.label, infoB.description.label);
        CPPUNIT_ASSERT_EQUAL(_data->descriptionB.vectorFieldType, infoB.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(_data->descriptionB.scale, infoB.description.scale);
        const string_vector& componentNames = infoB.description.componentNames;
        CPPUNIT_ASSERT_EQUAL(_data->descriptionB.numComponents, componentNames.size());
        for (size_t i = 0; i < _data->descriptionB.numComponents; ++i) {
            CPPUNIT_ASSERT_EQUAL(_data->descriptionB.componentNames[i], componentNames[i]);
        } // for
        CPPUNIT_ASSERT_EQUAL(_data->discretizationB.basisOrder, infoB.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(_data->discretizationB.quadOrder, infoB.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(_data->discretizationB.feSpace, infoB.fe.feSpace);
        CPPUNIT_ASSERT_EQUAL(_data->discretizationB.isBasisContinuous, infoB.fe.isBasisContinuous);
    } // Test getSubfieldInfo() for subfield B.

    CPPUNIT_ASSERT_THROW(_field->getSubfieldInfo("aabbccdd"), std::runtime_error);

    PYLITH_METHOD_END;
} /// testSubfieldAccessors


// ---------------------------------------------------------------------------------------------------------------------
// Test allocate().
void
pylith::topology::TestFieldMesh::testAllocate(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);

    _checkValues(*_field);

    PYLITH_METHOD_END;
} // testAllocate


// ---------------------------------------------------------------------------------------------------------------------
// Test zeroLocal().
void
pylith::topology::TestFieldMesh::testZeroLocal(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_field);
    _field->zeroLocal();

    _checkValues(*_field, 0.0);

    PYLITH_METHOD_END;
} // testZeroLocal


// ---------------------------------------------------------------------------------------------------------------------
// Test view().
void
pylith::topology::TestFieldMesh::testView(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_field);
    _field->view("Testing view");

    PYLITH_METHOD_END;
} // testView


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    const int cellDim = _data->cellDim;
    const int numCells = _data->numCells;
    const int numVertices = _data->numVertices;
    const int numCorners = _data->numCorners;
    const int spaceDim = _data->cellDim;

    PylithInt size = numVertices * spaceDim;
    scalar_array coordinates(size);
    for (PylithInt i = 0; i < size; ++i) {
        coordinates[i] = _data->coordinates[i];
    } // for

    size = numCells * numCorners;
    int_array cells(size);
    for (PylithInt i = 0; i < size; ++i) {
        cells[i] = _data->cells[i];
    } // for

    delete _mesh;_mesh = new Mesh;CPPUNIT_ASSERT(_mesh);
    pylith::meshio::MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners,
                                           cellDim);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(spaceDim);
    _mesh->setCoordSys(&cs);

    // Setup labels for constraints.
    PetscErrorCode err;
    for (PylithInt i = 0; i < _data->bcANumVertices; ++i) {
        err = DMSetLabelValue(_mesh->getDM(), _data->bcALabel, numCells+_data->bcAVertices[i], _data->bcALabelId);
        CPPUNIT_ASSERT(!err);
    } // for
    for (PylithInt i = 0; i < _data->bcBNumVertices; ++i) {
        err = DMSetLabelValue(_mesh->getDM(), _data->bcBLabel, numCells+_data->bcBVertices[i], _data->bcBLabelId);
        CPPUNIT_ASSERT(!err);
    } // for

    // Setup field
    delete _field;_field = new Field(*_mesh);
    _field->setLabel("solution");
    _field->subfieldAdd(_data->descriptionA, _data->discretizationA);
    _field->subfieldAdd(_data->descriptionB, _data->discretizationB);
    _field->subfieldsSetup();
    _field->createDiscretization();

    PetscDMLabel labelA = NULL, labelB = NULL;
    err = DMGetLabel(_field->getDM(), _data->bcALabel, &labelA);CPPUNIT_ASSERT(!err);
    err = DMGetLabel(_field->getDM(), _data->bcBLabel, &labelB);CPPUNIT_ASSERT(!err);
    const PetscInt numLabelValues = 1;
    PetscInt i_field = 0;
    err = DMAddBoundary(_field->getDM(), DM_BC_ESSENTIAL, "bcA", labelA, numLabelValues, &_data->bcALabelId, i_field,
                        _data->bcANumConstrainedDOF, _data->bcAConstrainedDOF, NULL, NULL, NULL, NULL);CPPUNIT_ASSERT(!err);
    i_field = 1;
    err = DMAddBoundary(_field->getDM(), DM_BC_ESSENTIAL, "bcB", labelB, numLabelValues, &_data->bcBLabelId, i_field,
                        _data->bcBNumConstrainedDOF, _data->bcBConstrainedDOF, NULL, NULL, NULL, NULL);CPPUNIT_ASSERT(!err);
    // Allocate field.
    _field->allocate();

    // Populate with values.
    PetscDM dmMesh = _mesh->getDM();CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PylithInt vStart = depthStratum.begin();
    const PylithInt vEnd = depthStratum.end();

    VecVisitorMesh fieldVisitor(*_field);
    const PylithInt fiberDim = _data->descriptionA.numComponents + _data->descriptionB.numComponents;
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for (PylithInt v = vStart, indexA = 0, indexB = 0; v < vEnd; ++v) {
        // Set values for field A
        const PylithInt offA = fieldVisitor.sectionOffset(v);
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
        for (size_t d = 0; d < _data->descriptionA.numComponents; ++d) {
            fieldArray[offA+d] = _data->subfieldAValues[indexA++];
        } // for
          // Set values for field B
        const PylithInt offB = offA + _data->descriptionA.numComponents;
        for (size_t d = 0; d < _data->descriptionB.numComponents; ++d) {
            fieldArray[offB+d] = _data->subfieldBValues[indexB++];
        } // for
    } // for

    PYLITH_METHOD_END;
} // _initialize


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_checkValues(const Field& field,
                                              const PylithReal scale) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    // Create array of values in field from subfields.
    const int numVertices = _data->numVertices;
    const int fiberDimA = _data->descriptionA.numComponents;
    const int fiberDimB = _data->descriptionB.numComponents;
    scalar_array valuesE(numVertices * (fiberDimA + fiberDimB));
    for (int iVertex = 0, index = 0, indexA = 0, indexB = 0; iVertex < numVertices; ++iVertex) {
        for (int d = 0; d < fiberDimA; ++d) {
            valuesE[index++] = _data->subfieldAValues[indexA++];
        } // for
        for (int d = 0; d < fiberDimB; ++d) {
            valuesE[index++] = _data->subfieldBValues[indexB++];
        } // for
    } // for

    PetscDM dmMesh = _mesh->getDM();CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PylithInt vStart = depthStratum.begin();
    const PylithInt vEnd = depthStratum.end();

    VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    const PylithInt fiberDim = fiberDimA + fiberDimB;
    const PylithReal tolerance = 1.0e-6;
    for (PylithInt v = vStart, index = 0; v < vEnd; ++v) {
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
        const PylithInt off = fieldVisitor.sectionOffset(v);

        for (PylithInt d = 0; d < fiberDim; ++d) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[index++]*scale, fieldArray[off+d], tolerance);
        } // for
    } // for

    PYLITH_METHOD_END;
} // _checkValues


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_checkValues(const PetscVec& vec,
                                              const PylithReal scale) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    // Create array of values in field from subfields.
    const int numVertices = _data->numVertices;
    const int fiberDimA = _data->descriptionA.numComponents;
    const int fiberDimB = _data->descriptionB.numComponents;
    scalar_array valuesE(numVertices * (fiberDimA + fiberDimB));
    for (int iVertex = 0, index = 0, indexA = 0, indexB = 0; iVertex < numVertices; ++iVertex) {
        for (int d = 0; d < fiberDimA; ++d) {
            valuesE[index++] = _data->subfieldAValues[indexA++];
        } // for
        for (int d = 0; d < fiberDimB; ++d) {
            valuesE[index++] = _data->subfieldBValues[indexB++];
        } // for
    } // for

    PetscErrorCode err;
    PylithInt size = 0;
    PylithScalar* vecArray = NULL;
    err = VecGetSize(vec, &size);CPPUNIT_ASSERT(!err);
    err = VecGetArray(vec, &vecArray);CPPUNIT_ASSERT(!err);

    const PylithInt sizeE = numVertices * (fiberDimA + fiberDimB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_EQUAL(sizeE, size);
    for (PylithInt i = 0; i < sizeE; ++i) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i]*scale, vecArray[i], tolerance);
    } // for
    err = VecRestoreArray(vec, &vecArray);CPPUNIT_ASSERT(!err);

    PYLITH_METHOD_END;
} // _checkValues


// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::topology::TestFieldMesh_Data::TestFieldMesh_Data(void) :
    cellDim(0),
    numVertices(0),
    numCells(0),
    numCorners(0),
    cells(NULL),
    coordinates(NULL),

    subfieldAValues(NULL),
    bcALabel(NULL),
    bcALabelId(0),
    bcANumConstrainedDOF(0),
    bcAConstrainedDOF(NULL),
    bcANumVertices(0),
    bcAVertices(NULL),

    subfieldBValues(NULL),
    bcBLabel(NULL),
    bcBLabelId(0),
    bcBNumConstrainedDOF(0),
    bcBConstrainedDOF(NULL),
    bcBNumVertices(0),
    bcBVertices(NULL) { // constructor
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::topology::TestFieldMesh_Data::~TestFieldMesh_Data(void) {}


// End of file
