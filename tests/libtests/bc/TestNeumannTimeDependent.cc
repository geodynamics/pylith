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

#include "TestNeumannTimeDependent.hh" // Implementation of class methods

#include "pylith/bc/NeumannTimeDependent.hh" // USES NeumannTimeDependent

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

const double pylith::bc::TestNeumannTimeDependent::FILL_VALUE = -999.0;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestNeumannTimeDependent::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _bc = new pylith::bc::NeumannTimeDependent();CPPUNIT_ASSERT(_bc);

    _data = NULL;
    _mesh = NULL;
    _solution = NULL;

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestNeumannTimeDependent::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _bc;_bc = NULL;

    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;
    delete _solution;_solution = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestNeumannTimeDependent::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    NeumannTimeDependent* bc = new NeumannTimeDependent();CPPUNIT_ASSERT(bc);
    delete bc;bc = NULL;

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
/// Test accessors (dbTimeHistory, useInitial, useRate, useTimeHistory).
void
pylith::bc::TestNeumannTimeDependent::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_bc);

    bool flag;

    // field()
    const std::string fieldName = "displacement";
    _bc->field(fieldName.c_str());
    CPPUNIT_ASSERT_EQUAL(fieldName, std::string(_bc->field()));

    // dbTimeHistory()
    spatialdata::spatialdb::TimeHistory db;
    CPPUNIT_ASSERT_EQUAL((void*)NULL, (void*)_bc->dbTimeHistory());
    _bc->dbTimeHistory(&db);
    CPPUNIT_ASSERT_EQUAL((void*)&db, (void*)_bc->dbTimeHistory());

    // useInitial()
    CPPUNIT_ASSERT_EQUAL(true, _bc->useInitial());
    flag = true;
    _bc->useInitial(flag);
    CPPUNIT_ASSERT_EQUAL(flag, _bc->useInitial());
    flag = false;
    _bc->useInitial(flag);
    CPPUNIT_ASSERT_EQUAL(flag, _bc->useInitial());

    // useRate()
    CPPUNIT_ASSERT_EQUAL(false, _bc->useRate());
    flag = true;
    _bc->useRate(flag);
    CPPUNIT_ASSERT_EQUAL(flag, _bc->useRate());
    flag = false;
    _bc->useRate(flag);
    CPPUNIT_ASSERT_EQUAL(flag, _bc->useRate());

    // useTimeHistory()
    CPPUNIT_ASSERT_EQUAL(false, _bc->useTimeHistory());
    flag = true;
    _bc->useTimeHistory(flag);
    CPPUNIT_ASSERT_EQUAL(flag, _bc->useTimeHistory());
    flag = false;
    _bc->useTimeHistory(flag);
    CPPUNIT_ASSERT_EQUAL(flag, _bc->useTimeHistory());

    PYLITH_METHOD_END;
} // testAccessors


// ----------------------------------------------------------------------
void
pylith::bc::TestNeumannTimeDependent::testAuxFieldDiscretization(void) {
    PYLITH_METHOD_BEGIN;

    const topology::FieldBase::Discretization infoDefault = pylith::topology::Field::Discretization(1, 1, true, pylith::topology::FieldBase::POLYNOMIAL_SPACE);
    const topology::FieldBase::Discretization infoA = pylith::topology::Field::Discretization(1, 2, false, pylith::topology::FieldBase::POLYNOMIAL_SPACE);
    const topology::FieldBase::Discretization infoB = pylith::topology::Field::Discretization(2, 2, true, pylith::topology::FieldBase::POINT_SPACE);

    CPPUNIT_ASSERT(_bc);
    _bc->auxSubfieldDiscretization("A", infoA.basisOrder, infoA.quadOrder, infoA.isBasisContinuous, infoA.feSpace);
    _bc->auxSubfieldDiscretization("B", infoB.basisOrder, infoB.quadOrder, infoB.isBasisContinuous, infoB.feSpace);

    CPPUNIT_ASSERT(_bc->_auxiliaryFactory());
    { // A
        const topology::FieldBase::Discretization& test = _bc->_auxiliaryFactory()->getSubfieldDiscretization("A");
        CPPUNIT_ASSERT_EQUAL(infoA.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoA.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoA.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoA.feSpace, test.feSpace);
    } // A

    { // B
        const topology::FieldBase::Discretization& test = _bc->_auxiliaryFactory()->getSubfieldDiscretization("B");
        CPPUNIT_ASSERT_EQUAL(infoB.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoB.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoB.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoB.feSpace, test.feSpace);
    } // B

    { // C (default)
        const topology::FieldBase::Discretization& test = _bc->_auxiliaryFactory()->getSubfieldDiscretization("C");
        CPPUNIT_ASSERT_EQUAL(infoDefault.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoDefault.feSpace, test.feSpace);
    } // C (default)

    { // default
        const topology::FieldBase::Discretization& test = _bc->_auxiliaryFactory()->getSubfieldDiscretization("default");
        CPPUNIT_ASSERT_EQUAL(infoDefault.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoDefault.feSpace, test.feSpace);
    } // default

    PYLITH_METHOD_END;
} // testAuxFieldDiscretization


// ----------------------------------------------------------------------
// Test auxFieldDB().
void
pylith::bc::TestNeumannTimeDependent::testAuxFieldDB(void) {
    PYLITH_METHOD_BEGIN;

    const std::string label = "test db";
    spatialdata::spatialdb::UserFunctionDB db;
    db.setLabel(label.c_str());

    CPPUNIT_ASSERT(_bc);
    _bc->auxFieldDB(&db);

    CPPUNIT_ASSERT(_bc->_auxiliaryFactory());
    CPPUNIT_ASSERT(_bc->_auxiliaryFactory()->queryDB());

    CPPUNIT_ASSERT_EQUAL(label, std::string(_bc->_auxiliaryFactory()->queryDB()->getLabel()));

    PYLITH_METHOD_END;
} // testAuxFieldDB


// ----------------------------------------------------------------------
// Test normalizer().
void
pylith::bc::TestNeumannTimeDependent::testNormalizer(void) {
    PYLITH_METHOD_BEGIN;

    spatialdata::units::Nondimensional normalizer;
    const double scale = 5.0;
    normalizer.setLengthScale(scale);

    CPPUNIT_ASSERT(_bc);
    _bc->normalizer(normalizer);
    CPPUNIT_ASSERT_EQUAL(scale, _bc->_normalizer->getLengthScale());

    PYLITH_METHOD_END;
} // testNormalizer


// ----------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::bc::TestNeumannTimeDependent::testVerifyConfiguration(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    _setupSolutionField();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);

    // Verify should pass.
    _bc->verifyConfiguration(*_solution);

    // Check for failure with field not in solution.
    _bc->field("dslfjadsf");
    CPPUNIT_ASSERT_THROW(_bc->verifyConfiguration(*_solution), std::runtime_error);

    PYLITH_METHOD_END;
} // testVerifyConfiguration


// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestNeumannTimeDependent::testInitialize(void) {
    PYLITH_METHOD_BEGIN;

    // Call initialize()
    _initialize(); // includes setting up auxField
    _setupSolutionField();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);

    _bc->initialize(*_solution);

    // _bc->auxField().view("AUXILIARY FIELD"); // DEBUGGING

    // Verify auxiliary field.
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_mesh);
    const pylith::topology::Field* auxField = _bc->auxField();CPPUNIT_ASSERT(auxField);
    CPPUNIT_ASSERT_EQUAL(std::string("auxiliary subfields"), std::string(auxField->getLabel()));
    CPPUNIT_ASSERT_EQUAL(_mesh->dimension(), auxField->getSpaceDim());

    PylithReal norm = 0.0;
    PylithReal t = _data->t;
    const PetscDM dm = auxField->dmMesh();CPPUNIT_ASSERT(dm);
    pylith::topology::FieldQuery query(*auxField);
    query.initializeWithDefaultQueryFns();
    CPPUNIT_ASSERT(_data->normalizer);
    query.openDB(_data->auxDB, _data->normalizer->getLengthScale());
    PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, query.functions(), (void**)query.contextPtrs(), auxField->localVector(), &norm);CPPUNIT_ASSERT(!err);
    query.closeDB(_data->auxDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testInitialize


// ----------------------------------------------------------------------
// Test preStep().
void
pylith::bc::TestNeumannTimeDependent::testPrestep(void) {
    PYLITH_METHOD_BEGIN;

    if (!_data->useTimeHistory) {PYLITH_METHOD_END;}

    _initialize();
    _setupSolutionField();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);
    _bc->initialize(*_solution);

    const pylith::topology::Field* auxField = _bc->auxField();CPPUNIT_ASSERT(auxField);
    CPPUNIT_ASSERT(auxField->hasSubfield("time_history_value"));

    pylith::topology::Field valueField(*_mesh);
    valueField.copySubfield(*auxField, "time_history_value");
    CPPUNIT_ASSERT(valueField.sectionSize() > 0);
    CPPUNIT_ASSERT_EQUAL(std::string("time_history_value"), std::string(valueField.getLabel()));

    PylithReal norm = 0.0;
    PylithReal t = _data->t;
    const PetscDM dm = auxField->dmMesh();CPPUNIT_ASSERT(dm);
    pylith::topology::FieldQuery query(valueField);
    query.initializeWithDefaultQueryFns();
    CPPUNIT_ASSERT(_data->normalizer);
    query.openDB(_data->auxDB, _data->normalizer->getLengthScale());
    PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, query.functions(), (void**)query.contextPtrs(), valueField.localVector(), &norm);CPPUNIT_ASSERT(!err);
    query.closeDB(_data->auxDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testPrestep


// ----------------------------------------------------------------------
// Test setSolution().
void
pylith::bc::TestNeumannTimeDependent::testComputeRHSResidual(void) {
    PYLITH_METHOD_BEGIN;

#if 1
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad not implemented.", false);
#else
    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);
    _solution->createDiscretization();
    _bc->initialize(*_solution);

    // Initialize solution field.
    _solution->allocate();
    _solution->createScatter(_solution->mesh(), "global");

    // Set solution field.
    CPPUNIT_ASSERT(_data);
    _bc->setSolution(*_solution, _data->t);

    // Verify number and DOF of constraints in solution field.
    int iConstraint = 0;
    PetscErrorCode err = 0;
    for (PetscInt v = vStart; v < vEnd; ++v) {
        PetscInt dof, cdof, fdof, fcdof;

        err = PetscSectionGetDof(fieldSection, v, &dof);PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetConstraintDof(fieldSection, v, &cdof);PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetFieldDof(fieldSection, v, 0, &fdof);PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetFieldConstraintDof(fieldSection, v, 0, &fcdof);PYLITH_CHECK_ERROR(err);
        if (v != _data->constrainedPoints[iConstraint] + offset) {
            CPPUNIT_ASSERT_EQUAL(_data->numDOF, dof);
            CPPUNIT_ASSERT_EQUAL(0, cdof);
            CPPUNIT_ASSERT_EQUAL(_data->numDOF, fdof);
            CPPUNIT_ASSERT_EQUAL(0, fcdof);
        } else {
            CPPUNIT_ASSERT_EQUAL(_data->numDOF, dof);
            CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, cdof);
            CPPUNIT_ASSERT_EQUAL(_data->numDOF, fdof);
            CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, fcdof);
            ++iConstraint;
        } // if/else
    } // for

    // Verify values in solution field.
    PylithReal norm = 0.0;
    PylithReal t = _data->t;
    const PetscDM dmSoln = _solution->dmMesh();CPPUNIT_ASSERT(dmSoln);
    pylith::topology::FieldQuery* query = _db->_auxiliaryFieldsQuery;
    query->openDB(queryDB, _data->lengthScale);

    PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, query->functions(), (void**)query->contextPtrs(), _solution->localVector(), &norm);CPPUNIT_ASSERT(!err);
    query->closeDB(queryDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
#endif

    PYLITH_METHOD_END;
} // testComputeRHSResidual


// ----------------------------------------------------------------------
// Test _auxiliaryFieldSetup().
void
pylith::bc::TestNeumannTimeDependent::testAuxFieldSetup(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    _setupSolutionField();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);
    const PylithReal timeScale = _data->normalizer->getTimeScale();

    delete _bc->_boundaryMesh;_bc->_boundaryMesh = new pylith::topology::Mesh(_solution->mesh(), _data->bcLabel);
    CPPUNIT_ASSERT(_bc->_boundaryMesh);

    delete _bc->_auxiliaryField;_bc->_auxiliaryField = new pylith::topology::Field(*_bc->_boundaryMesh);CPPUNIT_ASSERT(_bc->_auxiliaryField);
    _bc->_auxiliaryFieldSetup(*_solution);

    CPPUNIT_ASSERT(_mesh->getCoordSys());
    const size_t spaceDim = _mesh->getCoordSys()->getSpaceDim();
    const pylith::topology::Field::VectorFieldEnum vectorFieldType = _data->vectorFieldType;
    const size_t numComponents = (vectorFieldType == pylith::topology::Field::VECTOR) ? spaceDim : 1;

    // Check discretizations
    int ifield = 0;
    if (_data->useInitial) {
        const char* label = "initial_amplitude";
        CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(_data->auxSubfields[ifield]));
        const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];

        const pylith::topology::Field::SubfieldInfo& info = _bc->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(numComponents, info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(vectorFieldType, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(_data->scale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(discretization.basisOrder, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(discretization.quadOrder, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(discretization.isBasisContinuous, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(discretization.feSpace, info.fe.feSpace);
        ++ifield;
    } // if

    if (_data->useRate) {
        { // amplitude
            const char* label = "rate_amplitude";
            CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(_data->auxSubfields[ifield]));
            const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];

            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxiliaryField->subfieldInfo(label);
            CPPUNIT_ASSERT_EQUAL(numComponents, info.description.numComponents);
            CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
            CPPUNIT_ASSERT_EQUAL(vectorFieldType, info.description.vectorFieldType);
            CPPUNIT_ASSERT_EQUAL(_data->scale / timeScale, info.description.scale);
            CPPUNIT_ASSERT_EQUAL(discretization.basisOrder, info.fe.basisOrder);
            CPPUNIT_ASSERT_EQUAL(discretization.quadOrder, info.fe.quadOrder);
            CPPUNIT_ASSERT_EQUAL(discretization.isBasisContinuous, info.fe.isBasisContinuous);
            CPPUNIT_ASSERT_EQUAL(discretization.feSpace, info.fe.feSpace);
            ++ifield;
        } // amplitude

        { // start time
            const char* label = "rate_start_time";
            CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(_data->auxSubfields[ifield]));
            const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];

            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxiliaryField->subfieldInfo(label);
            CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
            CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
            CPPUNIT_ASSERT_EQUAL(timeScale, info.description.scale);
            CPPUNIT_ASSERT_EQUAL(discretization.basisOrder, info.fe.basisOrder);
            CPPUNIT_ASSERT_EQUAL(discretization.quadOrder, info.fe.quadOrder);
            CPPUNIT_ASSERT_EQUAL(discretization.isBasisContinuous, info.fe.isBasisContinuous);
            CPPUNIT_ASSERT_EQUAL(discretization.feSpace, info.fe.feSpace);
            ++ifield;
        } // start time
    } // if

    if (_data->useTimeHistory) {
        { // amplitude
            const char* label = "time_history_amplitude";
            CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(_data->auxSubfields[ifield]));
            const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];

            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxiliaryField->subfieldInfo(label);
            CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
            CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
            CPPUNIT_ASSERT_EQUAL(_data->scale, info.description.scale);
            CPPUNIT_ASSERT_EQUAL(discretization.basisOrder, info.fe.basisOrder);
            CPPUNIT_ASSERT_EQUAL(discretization.quadOrder, info.fe.quadOrder);
            CPPUNIT_ASSERT_EQUAL(discretization.isBasisContinuous, info.fe.isBasisContinuous);
            CPPUNIT_ASSERT_EQUAL(discretization.feSpace, info.fe.feSpace);
            ++ifield;
        } // amplitude

        { // start time
            const char* label = "time_history_start_time";
            CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(_data->auxSubfields[ifield]));
            const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];

            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxiliaryField->subfieldInfo(label);
            CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
            CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
            CPPUNIT_ASSERT_EQUAL(timeScale, info.description.scale);
            CPPUNIT_ASSERT_EQUAL(discretization.basisOrder, info.fe.basisOrder);
            CPPUNIT_ASSERT_EQUAL(discretization.quadOrder, info.fe.quadOrder);
            CPPUNIT_ASSERT_EQUAL(discretization.isBasisContinuous, info.fe.isBasisContinuous);
            CPPUNIT_ASSERT_EQUAL(discretization.feSpace, info.fe.feSpace);
            ++ifield;
        } // start time

        { // time_history_amplitude
            const char* label = "time_history_value";
            CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(_data->auxSubfields[ifield]));
            const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];

            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxiliaryField->subfieldInfo(label);
            CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
            CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
            CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
            CPPUNIT_ASSERT_EQUAL(discretization.basisOrder, info.fe.basisOrder);
            CPPUNIT_ASSERT_EQUAL(discretization.quadOrder, info.fe.quadOrder);
            CPPUNIT_ASSERT_EQUAL(discretization.isBasisContinuous, info.fe.isBasisContinuous);
            CPPUNIT_ASSERT_EQUAL(discretization.feSpace, info.fe.feSpace);
            ++ifield;
        } // time_history_amplitude

    } // if

    PYLITH_METHOD_END;
} // testAuxFieldSetup


// ----------------------------------------------------------------------
void
pylith::bc::TestNeumannTimeDependent::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_data);

    delete _mesh;_mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    pylith::meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(_data->meshFilename);
    iohandler.filename(_data->meshFilename);
    iohandler.read(_mesh);
    _mesh->setCoordSys(_data->cs);
    CPPUNIT_ASSERT(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    _bc->setLabel(_data->bcLabel);
    _bc->field(_data->field);
    _bc->auxFieldDB(_data->auxDB);
    _bc->normalizer(*_data->normalizer);
    for (int ifield = 0; ifield < _data->numAuxSubfields; ++ifield) {
        const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];
        const char* name = _data->auxSubfields[ifield];
        _bc->auxSubfieldDiscretization(name, discretization.basisOrder, discretization.quadOrder, discretization.isBasisContinuous, discretization.feSpace);
    } // for

    _bc->useInitial(_data->useInitial);
    _bc->useRate(_data->useRate);
    _bc->useTimeHistory(_data->useTimeHistory);
    if (_data->useTimeHistory) {
        CPPUNIT_ASSERT(_data->thFilename);
        spatialdata::spatialdb::TimeHistory th;
        th.filename(_data->thFilename);
        _bc->dbTimeHistory(&th);
    } // if

    PYLITH_METHOD_END;
} // _initialize


// ----------------------------------------------------------------------
void
pylith::bc::TestNeumannTimeDependent::_setupSolutionField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);

    delete _solution;_solution = new pylith::topology::Field(*_mesh);
    _solution->setLabel("Solution (displacement, velocity, pressure)");
    pylith::problems::SolutionFactory factory(*_solution, *_data->normalizer);
    factory.displacement(_data->solnDiscretizations[0]);
    factory.velocity(_data->solnDiscretizations[1]);
    factory.fluidPressure(_data->solnDiscretizations[2]);
    _solution->subfieldsSetup();

    PYLITH_METHOD_END;
} // setupSolutionField


// ----------------------------------------------------------------------
// Constructor
pylith::bc::TestNeumannTimeDependent_Data::TestNeumannTimeDependent_Data(void) :
    meshFilename(NULL),
    bcLabel(NULL),
    normalizer(new spatialdata::units::Nondimensional),
    field(NULL),
    vectorFieldType(pylith::topology::Field::OTHER),
    useInitial(false),
    useRate(false),
    useTimeHistory(false),
    thFilename(NULL),
    numAuxSubfields(0),
    auxSubfields(NULL),
    auxDiscretizations(NULL),
    auxDB(new spatialdata::spatialdb::UserFunctionDB),
    t(0.0),
    solnNumSubfields(0),
    solnDiscretizations(NULL),
    solnDB(new spatialdata::spatialdb::UserFunctionDB) {}


// ----------------------------------------------------------------------
// Destructor
pylith::bc::TestNeumannTimeDependent_Data::~TestNeumannTimeDependent_Data(void) {
    delete cs;cs = NULL;
    delete normalizer;normalizer = NULL;
    delete auxDB;auxDB = NULL;
    delete solnDB;solnDB = NULL;
} // destructor


// End of file
