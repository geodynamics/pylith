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
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory
#include "pylith/utils/constdefs.h" // USES PYLITH_MAXFLOAT

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

const double pylith::bc::TestDirichletTimeDependent::FILL_VALUE = -999.0;

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
    _bc->auxSubfieldDiscretization("A", infoA.basisOrder, infoA.quadOrder, infoA.isBasisContinuous, infoA.feSpace);
    _bc->auxSubfieldDiscretization("B", infoB.basisOrder, infoB.quadOrder, infoB.isBasisContinuous, infoB.feSpace);

    CPPUNIT_ASSERT(_bc->_auxFactory());
    { // A
        const topology::FieldBase::Discretization& test = _bc->_auxFactory()->subfieldDiscretization("A");
        CPPUNIT_ASSERT_EQUAL(infoA.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoA.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoA.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoA.feSpace, test.feSpace);
    } // A

    { // B
        const topology::FieldBase::Discretization& test = _bc->_auxFactory()->subfieldDiscretization("B");
        CPPUNIT_ASSERT_EQUAL(infoB.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoB.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoB.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoB.feSpace, test.feSpace);
    } // B

    { // C (default)
        const topology::FieldBase::Discretization& test = _bc->_auxFactory()->subfieldDiscretization("C");
        CPPUNIT_ASSERT_EQUAL(infoDefault.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoDefault.feSpace, test.feSpace);
    } // C (default)

    { // default
        const topology::FieldBase::Discretization& test = _bc->_auxFactory()->subfieldDiscretization("default");
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
pylith::bc::TestDirichletTimeDependent::testAuxFieldDB(void)
{ // testAuxFieldDB
    PYLITH_METHOD_BEGIN;

    const std::string label = "test db";
    spatialdata::spatialdb::UserFunctionDB db;
    db.label(label.c_str());

    CPPUNIT_ASSERT(_bc);
    _bc->auxFieldDB(&db);

    CPPUNIT_ASSERT(_bc->_auxFactory());
    CPPUNIT_ASSERT(_bc->_auxFactory()->queryDB());

    CPPUNIT_ASSERT_EQUAL(label, std::string(_bc->_auxFactory()->queryDB()->label()));

    PYLITH_METHOD_END;
} // testAuxFieldDB


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
    _setupSolutionField();

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

    // Call initialize()
    _initialize(); // includes setting up auxField
    _setupSolutionField();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);

    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(_solution->dmMesh(), &prob); CPPUNIT_ASSERT(!err);
    int numBCBefore = 0;
    err = PetscDSGetNumBoundary(prob, &numBCBefore); CPPUNIT_ASSERT(!err);

    _bc->initialize(*_solution);

    //_bc->auxField().view("AUXILIARY FIELD"); // DEBUGGING

    // Verify auxiliary field.
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_mesh);
    const pylith::topology::Field& auxField = _bc->auxField();
    CPPUNIT_ASSERT_EQUAL(std::string("auxiliary subfields"), std::string(auxField.label()));
    CPPUNIT_ASSERT_EQUAL(_mesh->dimension(), auxField.spaceDim());

    PylithReal norm = 0.0;
    PylithReal t = _data->t;
    const PetscDM dm = auxField.dmMesh(); CPPUNIT_ASSERT(dm);
    pylith::topology::FieldQuery query(auxField);
    query.initializeWithDefaultQueryFns();
    CPPUNIT_ASSERT(_data->normalizer);
    query.openDB(_data->auxDB, _data->normalizer->lengthScale());
    err = DMPlexComputeL2DiffLocal(dm, t, query.functions(), (void**)query.contextPtrs(), auxField.localVector(), &norm); CPPUNIT_ASSERT(!err);
    query.closeDB(_data->auxDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

    // Verify boundary condition was added to DS.
    int numBCAfter = 0;
    err = PetscDSGetNumBoundary(prob, &numBCAfter); CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(numBCAfter, 1+numBCBefore);

    PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test preStep().
void
pylith::bc::TestDirichletTimeDependent::testPrestep(void)
{ // testPrestep
    PYLITH_METHOD_BEGIN;

    if (!_data->useTimeHistory) {PYLITH_METHOD_END;}

    _initialize();
    _setupSolutionField();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);
    _bc->initialize(*_solution);

    const pylith::topology::Field& auxField = _bc->auxField();
    CPPUNIT_ASSERT(auxField.hasSubfield("time_history_value"));

    pylith::topology::Field valueField(*_mesh);
    valueField.copySubfield(auxField, "time_history_value");
    CPPUNIT_ASSERT(valueField.sectionSize() > 0);
    CPPUNIT_ASSERT_EQUAL(std::string("time_history_value"), std::string(valueField.label()));

    PylithReal norm = 0.0;
    PylithReal t = _data->t;
    const PetscDM dm = auxField.dmMesh(); CPPUNIT_ASSERT(dm);
    pylith::topology::FieldQuery query(valueField);
    query.initializeWithDefaultQueryFns();
    CPPUNIT_ASSERT(_data->normalizer);
    query.openDB(_data->auxDB, _data->normalizer->lengthScale());
    PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, query.functions(), (void**)query.contextPtrs(), valueField.localVector(), &norm); CPPUNIT_ASSERT(!err);
    query.closeDB(_data->auxDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testPrestep

// ----------------------------------------------------------------------
// Test setSolution().
void
pylith::bc::TestDirichletTimeDependent::testSetSolution(void)
{ // testSetSolution
    PYLITH_METHOD_BEGIN;

    _initialize();
    _setupSolutionField();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);
    _bc->initialize(*_solution);

    // Initialize solution field.
    _solution->allocate();
    PetscErrorCode err;
    err = VecSet(_solution->localVector(), FILL_VALUE); CPPUNIT_ASSERT(!err);

    // Set solution field.
    CPPUNIT_ASSERT(_data);
    _solution->mesh().view(":detail.txt:ascii_info_detail");
    _bc->setSolution(_solution, _data->t);

    //_solution->view("SOLUTION BC ONLY"); // DEBUGGING

    // Verify setting solution did not change unconstrained values.
    const PylithReal tolerance = 1.0e-6;
    _solution->createScatter(_solution->mesh(), "global");
    _solution->scatterLocalToContext("global", ADD_VALUES);
    PylithScalar value = 0.0;
    err = VecMax(_solution->scatterVector("global"), NULL, &value); CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(FILL_VALUE, value, tolerance);
    err = VecMin(_solution->scatterVector("global"), NULL, &value); CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(FILL_VALUE, value, tolerance);

    // Verify solution values match expected values.
    // Fill unconstrained values only.
    const PylithReal t = _data->t;
    const PetscDM dmSoln = _solution->dmMesh(); CPPUNIT_ASSERT(dmSoln);
    pylith::topology::FieldQuery query(*_solution);
    query.initializeWithDefaultQueryFns();
    CPPUNIT_ASSERT(_data->normalizer);
    query.openDB(_data->solnDB, _data->normalizer->lengthScale());
    err = DMProjectFunction(dmSoln, t, query.functions(), (void**)query.contextPtrs(), INSERT_VALUES, _solution->scatterVector("global")); CPPUNIT_ASSERT(!err);
    query.closeDB(_data->solnDB);
    _solution->scatterContextToLocal("global", INSERT_VALUES);

    //_solution->view("SOLUTION ALL"); // DEBUGGING

    PylithReal norm = 0.0;
    query.openDB(_data->solnDB, _data->normalizer->lengthScale());
    err = DMPlexComputeL2DiffLocal(dmSoln, t, query.functions(), (void**)query.contextPtrs(), _solution->localVector(), &norm); CPPUNIT_ASSERT(!err);
    query.closeDB(_data->solnDB);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testSetSolution

// ----------------------------------------------------------------------
// Test _auxFieldSetup().
void
pylith::bc::TestDirichletTimeDependent::testAuxFieldSetup(void)
{ // testAuxFieldSetup
    PYLITH_METHOD_BEGIN;

    _initialize();
    _setupSolutionField();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);
    const PylithReal timeScale = _data->normalizer->timeScale();

    delete _bc->_boundaryMesh; _bc->_boundaryMesh = new pylith::topology::Mesh(_solution->mesh(), _data->bcLabel);
    CPPUNIT_ASSERT(_bc->_boundaryMesh);

    delete _bc->_auxField; _bc->_auxField = new pylith::topology::Field(*_bc->_boundaryMesh); CPPUNIT_ASSERT(_bc->_auxField);
    _bc->_auxFieldSetup(*_solution);

    CPPUNIT_ASSERT(_mesh->coordsys());
    const size_t spaceDim = _mesh->coordsys()->spaceDim();
    const pylith::topology::Field::VectorFieldEnum vectorFieldType = _data->vectorFieldType;
    const size_t numComponents = (vectorFieldType == pylith::topology::Field::VECTOR) ? spaceDim : 1;

    // Check discretizations
    int ifield = 0;
    if (_data->useInitial) {
        const char* label = "initial_amplitude";
        CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(_data->auxSubfields[ifield]));
        const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];

        const pylith::topology::Field::SubfieldInfo& info = _bc->_auxField->subfieldInfo(label);
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

            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxField->subfieldInfo(label);
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

            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxField->subfieldInfo(label);
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

            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxField->subfieldInfo(label);
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

            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxField->subfieldInfo(label);
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

            const pylith::topology::Field::SubfieldInfo& info = _bc->_auxField->subfieldInfo(label);
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
pylith::bc::TestDirichletTimeDependent::_initialize(void)
{ // _initialize
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_data);

    delete _mesh; _mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    pylith::meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(_data->meshFilename);
    iohandler.filename(_data->meshFilename);
    iohandler.read(_mesh);
    _mesh->coordsys(_data->cs);
    CPPUNIT_ASSERT(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    _bc->label(_data->bcLabel);
    _bc->field(_data->field);
    _bc->auxFieldDB(_data->auxDB);
    _bc->constrainedDOF(_data->constrainedDOF, _data->numConstrainedDOF);
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
pylith::bc::TestDirichletTimeDependent::_setupSolutionField(void)
{ // setupSolutionField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);

    delete _solution; _solution = new pylith::topology::Field(*_mesh);
    pylith::problems::SolutionFactory factory(*_solution, *_data->normalizer);
    factory.displacement(_data->solnDiscretizations[0]);
    factory.velocity(_data->solnDiscretizations[1]);
    factory.fluidPressure(_data->solnDiscretizations[2]);
    _solution->subfieldsSetup();

    PYLITH_METHOD_END;
} // setupSolutionField

// ----------------------------------------------------------------------
// Constructor
pylith::bc::TestDirichletTimeDependent_Data::TestDirichletTimeDependent_Data(void) :
    meshFilename(NULL),
    bcLabel(NULL),
    cs(NULL),
    normalizer(new spatialdata::units::Nondimensional),
    field(NULL),
    vectorFieldType(pylith::topology::Field::OTHER),
    numConstrainedDOF(0),
    constrainedDOF(NULL),
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
    solnDB(new spatialdata::spatialdb::UserFunctionDB)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::bc::TestDirichletTimeDependent_Data::~TestDirichletTimeDependent_Data(void)
{ // destructor
    delete cs; cs = NULL;
    delete normalizer; normalizer = NULL;
    delete auxDB; auxDB = NULL;
    delete solnDB; solnDB = NULL;
} // destructor


// End of file
