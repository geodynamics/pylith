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

#include "TestDirichletTimeDependent.hh" // Implementation of class methods

#include "pylith/bc/DirichletTimeDependent.hh" // USES DirichletTimeDependent

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleGridDB.hh" // USES SimpleGridDB
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/UniformDB.hh" // USES UniformDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletTimeDependent::setUp(void)
{ // setUp
    PYLITH_METHOD_BEGIN;

    _bc = new pylith::bc::DirichletTimeDependent();CPPUNIT_ASSERT(_bc);

    _data = NULL;
    _mesh = NULL;
    _solution = NULL;

    PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletTimeDependent::tearDown(void)
{ // tearDown
    PYLITH_METHOD_BEGIN;

    delete _bc; _bc = NULL;

    delete _data; _data = NULL;
    delete _mesh; _mesh = NULL;
    delete _solution; _solution = NULL;

    PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestDirichletTimeDependent::testConstructor(void)
{ // testConstructor
    PYLITH_METHOD_BEGIN;

    _bc = new DirichletTimeDependent();CPPUNIT_ASSERT(_bc);

    PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
/// Test accessors (dbTimeHistory, useInitial, useRate, useTimeHistory).
void
pylith::bc::TestDirichletTimeDependent::testAccessors(void)
{ // testAccessors
    PYLITH_METHOD_BEGIN;

    _bc = new DirichletTimeDependent();CPPUNIT_ASSERT(_bc);
    bool flag;

    // field()
    const std::string fieldName = "displacement";
    _bc->field(fieldName.c_str());
    CPPUNIT_ASSERT_EQUAL(fieldName, std::string(_bc->field()));

    // constrainedDOF().
    const size_t numConstrained = 4;
    const int constrainedDOF[4] = { 0, 2, 3, 6 };
    _bc->constrainedDOF(constrainedDOF, numConstrained);
    const pylith::int_array& dof = _bc->constrainedDOF();
    CPPUNIT_ASSERT_EQUAL(numConstrained, dof.size());
    for (size_t i = 0; i < numConstrained; ++i) {
        CPPUNIT_ASSERT_EQUAL(constrainedDOF[i], dof[i]);
    } // for

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
pylith::bc::TestDirichletTimeDependent::testAuxFieldDiscretization(void)
{ // testAuxFieldDiscretization
    PYLITH_METHOD_BEGIN;

    const topology::FieldBase::Discretization infoDefault = {-1, -1, true, pylith::topology::FieldBase::POLYNOMIAL_SPACE};
    const topology::FieldBase::Discretization infoA = {1, 2, false, pylith::topology::FieldBase::POLYNOMIAL_SPACE};
    const topology::FieldBase::Discretization infoB = {2, 2, true, pylith::topology::FieldBase::POINT_SPACE};

    CPPUNIT_ASSERT(_bc);
    _bc->auxFieldDiscretization("A", infoA.basisOrder, infoA.quadOrder, infoA.isBasisContinuous, infoA.feSpace);
    _bc->auxFieldDiscretization("B", infoB.basisOrder, infoB.quadOrder, infoB.isBasisContinuous, infoB.feSpace);

    { // A
        const topology::FieldBase::Discretization& test = _bc->auxFieldDiscretization("A");
        CPPUNIT_ASSERT_EQUAL(infoA.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoA.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoA.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoA.feSpace, test.feSpace);
    } // A

    { // B
        const topology::FieldBase::Discretization& test = _bc->auxFieldDiscretization("B");
        CPPUNIT_ASSERT_EQUAL(infoB.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoB.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoB.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoB.feSpace, test.feSpace);
    } // B

    { // C (default)
        const topology::FieldBase::Discretization& test = _bc->auxFieldDiscretization("C");
        CPPUNIT_ASSERT_EQUAL(infoDefault.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoDefault.feSpace, test.feSpace);
    } // C (default)

    { // default
        const topology::FieldBase::Discretization& test = _bc->auxFieldDiscretization("default");
        CPPUNIT_ASSERT_EQUAL(infoDefault.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoDefault.feSpace, test.feSpace);
    } // default


    PYLITH_METHOD_END;
} // testAuxFieldDiscretization

// ----------------------------------------------------------------------
// Test auxFieldsDB().
void
pylith::bc::TestDirichletTimeDependent::testAuxFieldsDB(void)
{ // testAuxFieldsDB
    PYLITH_METHOD_BEGIN;

    const std::string label = "test db";
    spatialdata::spatialdb::SimpleDB db;
    db.label(label.c_str());

    CPPUNIT_ASSERT(_bc);
    _bc->auxFieldsDB(&db);

    CPPUNIT_ASSERT(_bc->_auxFieldsDB);
    CPPUNIT_ASSERT_EQUAL(label, std::string(_bc->_auxFieldsDB->label()));

    PYLITH_METHOD_END;
} // testAuxFieldsDB


// ----------------------------------------------------------------------
// Test normalizer().
void
pylith::bc::TestDirichletTimeDependent::testNormalizer(void)
{ // testNormalizer
    PYLITH_METHOD_BEGIN;

    spatialdata::units::Nondimensional normalizer;
    const double scale = 5.0;
    normalizer.lengthScale(scale);

    CPPUNIT_ASSERT(_bc);
    _bc->normalizer(normalizer);
    CPPUNIT_ASSERT_EQUAL(scale, _bc->_normalizer->lengthScale());

    PYLITH_METHOD_END;
} // testNormalizer

// ----------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::bc::TestDirichletTimeDependent::testVerifyConfiguration(void)
{ // testVerifyConfiguration
    PYLITH_METHOD_BEGIN;

    _initialize();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);

    // Verify should pass.
    _bc->verifyConfiguration(*_solution);

    // Check for failure with field not in solution.
    _bc->field("dslfjadsf");
    CPPUNIT_ASSERT_THROW(_bc->verifyConfiguration(*_solution), std::runtime_error);

    // Check for failure with constrained DOF not in solution.
    const size_t numConstrainedDOF = 1;
    const int constrainedDOF[1] = { 9999 };
    _bc->constrainedDOF(constrainedDOF, numConstrainedDOF);
    CPPUNIT_ASSERT_THROW(_bc->verifyConfiguration(*_solution), std::runtime_error);

    PYLITH_METHOD_END;
} // testVerifyConfiguration

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestDirichletTimeDependent::testInitialize(void)
{ // testInitialize
    PYLITH_METHOD_BEGIN;

#if 1
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad not implemented.", false);
#else
    _initialize();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);
    _bc->initialize(*_solution);

    // Verify auxiliary fields.
    for (int i = 0; i < _data->numAuxFields; ++i) {
        CPPUNIT_ASSERT(_bc->hasAuxField(_data->auxFields[i]));
    } // for
    CPPUNIT_ASSERT(!_bc->hasAuxField("sgkgflgkjf"));

    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_mesh);
    const pylith::topology::Field& auxFields = _bc->auxFields();
    CPPUNIT_ASSERT_EQUAL(std::string("auxiliary fields"), std::string(auxFields.label()));
    CPPUNIT_ASSERT_EQUAL(_mesh->dimension(), auxFields.spaceDim());

    PylithReal norm = 0.0;
    PylithReal t = _data->t;
    const PetscDM dm = auxFields.dmMesh(); CPPUNIT_ASSERT(dm);
    pylith::topology::FieldQuery* query = _db->_auxFieldsQuery;
    query->openDB(_auxFieldsDB, _data->lengthScale);

    PetscErrorCode err = DMComputeL2Diff(dm, t, query->functions(), (void**)query->contextPtrs(), auxFields.localVector(), &norm); CPPUNIT_ASSERT(!err);
    query->closeDB(_auxFieldsDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

    // Verify boundary was added to DM.
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad Verify boundary was added to DM.", false);
#endif

    PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test preStep().
void
pylith::bc::TestDirichletTimeDependent::testPrestep(void)
{ // testPrestep
    PYLITH_METHOD_BEGIN;

#if 1
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad not implemented.", false);
#else
    _initialize();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);
    _bc->initialize(*_solution);

    pylith::topology::Field valueField(*_mesh);
    _bc->getAuxField(&valueField, "time_history_value");
    CPPUNIT_ASSERT(valueField.sectionSize() > 0);
    CPPUNIT_ASSERT_EQUAL(std::string("value"), std::string(valueField.label()));

    PylithReal norm = 0.0;
    PylithReal t = _data->t;
    const PetscDM dm = auxFields.dmMesh(); CPPUNIT_ASSERT(dm);
    pylith::topology::FieldQuery* query = _db->_auxFieldsQuery;
    query->openDB(_auxFieldsDB, _data->lengthScale);

    PetscErrorCode err = DMComputeL2Diff(dm, t, query->functions(), (void**)query->contextPtrs(), valueField.localVector(), &norm); CPPUNIT_ASSERT(!err);
    query->closeDB(_auxFieldsDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
#endif

    PYLITH_METHOD_END;
} // testPrestep

// ----------------------------------------------------------------------
// Test setSolution().
void
pylith::bc::TestDirichletTimeDependent::testSetSolution(void)
{ // testSetSolution
    PYLITH_METHOD_BEGIN;

#if 1
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad not implemented.", false);
#else
    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);
    _bc->initialize(*_solution);

    // Initialize solution field.
    _solution->allocate();
    _solution->zeroLocal();
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
    const PetscDM dmSoln = _solution->dmMesh(); CPPUNIT_ASSERT(dmSoln);
    pylith::topology::FieldQuery* query = _db->_auxFieldsQuery;
    query->openDB(_auxFieldsDB, _data->lengthScale);

    PetscErrorCode err = DMComputeL2Diff(dm, t, query->functions(), (void**)query->contextPtrs(), _solution->localVector(), &norm); CPPUNIT_ASSERT(!err);
    query->closeDB(_auxFieldsDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
#endif

    PYLITH_METHOD_END;
} // testSetSolution

// ----------------------------------------------------------------------
// Test _auxFieldsSetup().
void
pylith::bc::TestDirichletTimeDependent::testAuxFieldsSetup(void)
{ // testAuxFieldsSetup
    PYLITH_METHOD_BEGIN;

#if 0

    _setupSolutionFields();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);
    _bc->initialize(*_solution);

    CPPUNIT_ASSERT(_mymaterial);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_mydata);
    CPPUNIT_ASSERT(_mydata->normalizer);
    const PylithReal densityScale = _mydata->normalizer->densityScale();
    const PylithReal lengthScale = _mydata->normalizer->lengthScale();
    const PylithReal timeScale = _mydata->normalizer->timeScale();
    const PylithReal pressureScale = _mydata->normalizer->pressureScale();
    const PylithReal forceScale = densityScale * lengthScale / (timeScale * timeScale);


    CPPUNIT_ASSERT(_bc->_boundaryMesh);
    CPPUNIT_ASSERT(_data);

    delete _bc->_auxFields; _bc->_auxFields = new pylith::topology::Field(*_bc->_boundaryMesh); CPPUNIT_ASSERT(_bc->_auxFields);
    delete _bc->_auxFieldsQuery; _bc->_auxFieldsQuery = new pylith::topology::FieldQuery(*_bc->_auxFields); CPPUNIT_ASSERT(_bc->_auxFieldsQuery);
    _bc->_auxFieldsSetup();

    CPPUNIT_ASSERT(_mesh->coordsys());
    const int spaceDim = _mesh->coordsys()->spaceDim();
    const int numComponents = (_data->isConstrainedFieldVector) ? spaceDim : 1;
    const pylith::topology::Field::VectorFieldEnum vectorFieldType = (_data->isConstrainedFieldVector) ? pylith::topology::Field::VECTOR : pylith::topology::Field::SCALAR;

    // Check discretizations
    if (_data->useInitial) {
        const char* label = "initial_amplitude";
        const pylith::topology::Field::SubfieldInfo& info = _bc->_auxFields->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(numComponents, info.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.metadata.label);
        CPPUNIT_ASSERT_EQUAL(vectorFieldType, info.metadata.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(_data->lengthScale, info.metadata.scale);
        CPPUNIT_ASSERT_EQUAL(false, info.metadata.dimsOkay);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // if

    if (_data->useRate) {
        { // amplitude
            const char* label = "rate_amplitude";
            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxFields->subfieldInfo(label);
            CPPUNIT_ASSERT_EQUAL(numComponents, info.numComponents);
            CPPUNIT_ASSERT_EQUAL(std::string(label), info.metadata.label);
            CPPUNIT_ASSERT_EQUAL(vectorFieldType, info.metadata.vectorFieldType);
            CPPUNIT_ASSERT_EQUAL(_data->lengthScale / _data->timeScale, info.metadata.scale);
            CPPUNIT_ASSERT_EQUAL(false, info.metadata.dimsOkay);
            CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
            CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
            CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
        } // amplitude

        { // start time
            const char* label = "rate_start_time";
            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxFields->subfieldInfo(label);
            CPPUNIT_ASSERT_EQUAL(1, info.numComponents);
            CPPUNIT_ASSERT_EQUAL(std::string(label), info.metadata.label);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.metadata.vectorFieldType);
            CPPUNIT_ASSERT_EQUAL(_data->timeScale, info.metadata.scale);
            CPPUNIT_ASSERT_EQUAL(false, info.metadata.dimsOkay);
            CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
            CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
            CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
        } // start time
    } // if

    if (_data->useTimeHistory) {
        { // amplitude
            const char* label = "time_history_amplitude";
            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxFields->subfieldInfo(label);
            CPPUNIT_ASSERT_EQUAL(1, info.numComponents);
            CPPUNIT_ASSERT_EQUAL(std::string(label), info.metadata.label);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.metadata.vectorFieldType);
            CPPUNIT_ASSERT_EQUAL(_data->lengthScale, info.metadata.scale);
            CPPUNIT_ASSERT_EQUAL(false, info.metadata.dimsOkay);
            CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
            CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
            CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
        } // amplitude

        { // start time
            const char* label = "time_history_start_time";
            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxFields->subfieldInfo(label);
            CPPUNIT_ASSERT_EQUAL(1, info.numComponents);
            CPPUNIT_ASSERT_EQUAL(std::string(label), info.metadata.label);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.metadata.vectorFieldType);
            CPPUNIT_ASSERT_EQUAL(_data->timeScale, info.metadata.scale);
            CPPUNIT_ASSERT_EQUAL(false, info.metadata.dimsOkay);
            CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
            CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
            CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
        } // start time

        { // time_history_amplitude
            const char* label = "time_history_value";
            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxFields->subfieldInfo(label);
            CPPUNIT_ASSERT_EQUAL(1, info.numComponents);
            CPPUNIT_ASSERT_EQUAL(std::string(label), info.metadata.label);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.metadata.vectorFieldType);
            CPPUNIT_ASSERT_EQUAL(1.0, info.metadata.scale);
            CPPUNIT_ASSERT_EQUAL(false, info.metadata.dimsOkay);
            CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
            CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
            CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
            CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
        } // time_history_amplitude

    } // if

    // Make sure DB query functions are set correctly.
    if (_data->useInitial) {
        CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _bc->_auxFieldsQuery->queryFn("initial_amplitude"));
    } // if
    if (_data->useRate) {
        CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _bc->_auxFieldsQuery->queryFn("rate_amplitude"));
        CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _bc->_auxFieldsQuery->queryFn("rate_start_time"));
    } // if
    if (_data->useTimeHistory) {
        CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _bc->_auxFieldsQuery->queryFn("time_history_amplitude"));
        CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _bc->_auxFieldsQuery->queryFn("time_history_start_time"));
        CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _bc->_auxFieldsQuery->queryFn("time_history_value"));
    } // if
#endif
    PYLITH_METHOD_END;
} // testAuxFieldsSetup


// ----------------------------------------------------------------------
void
pylith::bc::TestDirichletTimeDependent::_initialize(void)
{ // _initialize
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_data);

    delete _mesh; _mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    pylith::meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(_data->meshFilename);
    iohandler.filename(_data->meshFilename);
    iohandler.read(_mesh);

    // Setup coordinates.
    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_mesh->dimension());
    cs.initialize();
    _mesh->coordsys(&cs);
    CPPUNIT_ASSERT(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

#if 0
    spatialdata::spatialdb::SimpleDB auxFieldsDB("TestDirichletTimeDependent auxFields");
    spatialdata::spatialdb::SimpleIOAscii dbIO;
    dbIO.filename(_data->auxDBFilename);
    auxFieldsDB.ioHandler(&dbIO);
    auxFieldsDB.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);
#endif

    _bc->label(_data->bcLabel);
    //bc->auxFieldsDB(&auxFieldsDB);
    _bc->constrainedDOF(_data->constrainedDOF, _data->numConstrainedDOF);
    _bc->normalizer(*_data->normalizer);

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
pylith::bc::TestDirichletTimeDependent::_setupSolutionField(void)
{ // setupSolutionField
    PYLITH_METHOD_BEGIN;

    spatialdata::spatialdb::SimpleGridDB solutionDB;
    solutionDB.filename(_data->solnDBFilename);
    solutionDB.label("solution database");
    solutionDB.queryType(spatialdata::spatialdb::SimpleGridDB::LINEAR);

    pylith::problems::SolutionFactory factory(*_solution, *_data->normalizer);
    factory.displacement(_data->solnDiscretizations[0]);
    factory.velocity(_data->solnDiscretizations[1]);
    factory.fluidPressure(_data->solnDiscretizations[2]);
    _solution->subfieldsSetup();
    _solution->allocate();
    factory.setValues(&solutionDB);

    PYLITH_METHOD_END;
} // setupSolutionField

// ----------------------------------------------------------------------
// Constructor
pylith::bc::TestDirichletTimeDependent_Data::TestDirichletTimeDependent_Data(void) :
    meshFilename(NULL),
    bcLabel(NULL),
    normalizer(new spatialdata::units::Nondimensional),
    field(NULL),
    vectorFieldType(pylith::topology::Field::OTHER),
    numConstrainedDOF(0),
    constrainedDOF(NULL),
    useInitial(false),
    useRate(false),
    useTimeHistory(false),
    thFilename(NULL),
    t(0.0),
    solnNumFields(0),
    solnDiscretizations(NULL),
    solnDBFilename(NULL),
    numAuxFields(0),
    auxFields(NULL),
    auxDiscretizations(NULL),
    auxDBFilename(NULL)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::bc::TestDirichletTimeDependent_Data::~TestDirichletTimeDependent_Data(void)
{ // destructor
    delete normalizer; normalizer = NULL;
} // destructor


// End of file
