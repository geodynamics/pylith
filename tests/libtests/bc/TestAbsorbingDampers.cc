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

#include "TestAbsorbingDampers.hh" // Implementation of class methods

#include "pylith/bc/AbsorbingDampers.hh" // USES AbsorbingDampers

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

const double pylith::bc::TestAbsorbingDampers::FILL_VALUE = -999.0;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestAbsorbingDampers::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _bc = new pylith::bc::AbsorbingDampers();CPPUNIT_ASSERT(_bc);

    _data = NULL;
    _mesh = NULL;
    _solution = NULL;

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestAbsorbingDampers::tearDown(void) {
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
pylith::bc::TestAbsorbingDampers::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    AbsorbingDampers* bc = new AbsorbingDampers();CPPUNIT_ASSERT(bc);
    delete bc;bc = NULL;

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
/// Test accessors (dbTimeHistory, useInitial, useRate, useTimeHistory).
void
pylith::bc::TestAbsorbingDampers::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_bc);

    // field()
    const std::string fieldName = "displacement";
    _bc->field(fieldName.c_str());
    CPPUNIT_ASSERT_EQUAL(fieldName, std::string(_bc->field()));

    PYLITH_METHOD_END;
} // testAccessors


// ----------------------------------------------------------------------
void
pylith::bc::TestAbsorbingDampers::testAuxFieldDiscretization(void) {
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
pylith::bc::TestAbsorbingDampers::testAuxFieldDB(void) {
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
pylith::bc::TestAbsorbingDampers::testNormalizer(void) {
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
pylith::bc::TestAbsorbingDampers::testVerifyConfiguration(void) {
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

    // Check for failure with field not in solution.
    _bc->field("fluid_pressure");
    CPPUNIT_ASSERT_THROW(_bc->verifyConfiguration(*_solution), std::runtime_error);

    PYLITH_METHOD_END;
} // testVerifyConfiguration


// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestAbsorbingDampers::testInitialize(void) {
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
    CPPUNIT_ASSERT_EQUAL(_mesh->getDimension(), auxField->getSpaceDim());

    PylithReal norm = 0.0;
    PylithReal t = _data->t;
    const PetscDM dm = auxField->getDM();CPPUNIT_ASSERT(dm);
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
// Test computeRHSResidual().
void
pylith::bc::TestAbsorbingDampers::testComputeRHSResidual(void) {
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
    const PetscDM dmSoln = _solution->getDM();CPPUNIT_ASSERT(dmSoln);
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
pylith::bc::TestAbsorbingDampers::testAuxFieldSetup(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    _setupSolutionField();

    CPPUNIT_ASSERT(_bc);
    CPPUNIT_ASSERT(_solution);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);

    delete _bc->_boundaryMesh;_bc->_boundaryMesh = new pylith::topology::Mesh(_solution->mesh(), _data->bcLabel);
    CPPUNIT_ASSERT(_bc->_boundaryMesh);

    delete _bc->_auxiliaryField;_bc->_auxiliaryField = new pylith::topology::Field(*_bc->_boundaryMesh);CPPUNIT_ASSERT(_bc->_auxiliaryField);
    _bc->_auxiliaryFieldSetup(*_solution);

    // Check discretizations
    int ifield = 0;
    { // density
        const char* label = "density";
        CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(_data->auxSubfields[ifield]));
        const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];

        const pylith::topology::Field::SubfieldInfo& info = _bc->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(_data->normalizer->getDensityScale(), info.description.scale);
        CPPUNIT_ASSERT_EQUAL(discretization.basisOrder, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(discretization.quadOrder, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(discretization.isBasisContinuous, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(discretization.feSpace, info.fe.feSpace);
        ++ifield;
    } // density

    const PylithScalar velocityScale = _data->normalizer->getLengthScale() / _data->normalizer->getTimeScale();

    { // vp
        const char* label = "vp";
        CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(_data->auxSubfields[ifield]));
        const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];

        const pylith::topology::Field::SubfieldInfo& info = _bc->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(velocityScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(discretization.basisOrder, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(discretization.quadOrder, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(discretization.isBasisContinuous, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(discretization.feSpace, info.fe.feSpace);
        ++ifield;
    } // vp

    { // vs
        const char* label = "vs";
        CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(_data->auxSubfields[ifield]));
        const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];

        const pylith::topology::Field::SubfieldInfo& info = _bc->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(velocityScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(discretization.basisOrder, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(discretization.quadOrder, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(discretization.isBasisContinuous, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(discretization.feSpace, info.fe.feSpace);
        ++ifield;
    } // vs

    PYLITH_METHOD_END;
} // testAuxFieldSetup


// ----------------------------------------------------------------------
void
pylith::bc::TestAbsorbingDampers::_initialize(void) {
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

    PYLITH_METHOD_END;
} // _initialize


// ----------------------------------------------------------------------
void
pylith::bc::TestAbsorbingDampers::_setupSolutionField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);

    delete _solution;_solution = new pylith::topology::Field(*_mesh);
    pylith::problems::SolutionFactory factory(*_solution, *_data->normalizer);
    factory.displacement(_data->solnDiscretizations[0]);
    factory.velocity(_data->solnDiscretizations[1]);
    factory.fluidPressure(_data->solnDiscretizations[2]);
    _solution->subfieldsSetup();

    PYLITH_METHOD_END;
} // setupSolutionField


// ----------------------------------------------------------------------
// Constructor
pylith::bc::TestAbsorbingDampers_Data::TestAbsorbingDampers_Data(void) :
    meshFilename(NULL),
    bcLabel(NULL),
    normalizer(new spatialdata::units::Nondimensional),
    field(NULL),
    vectorFieldType(pylith::topology::Field::OTHER),
    numAuxSubfields(0),
    auxSubfields(NULL),
    auxDiscretizations(NULL),
    auxDB(new spatialdata::spatialdb::UserFunctionDB),
    t(0.0),
    solnNumSubfields(0),
    solnDiscretizations(NULL),
    solnDB(new spatialdata::spatialdb::UserFunctionDB) {} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::bc::TestAbsorbingDampers_Data::~TestAbsorbingDampers_Data(void) {
    delete cs;cs = NULL;
    delete normalizer;normalizer = NULL;
    delete auxDB;auxDB = NULL;
    delete solnDB;solnDB = NULL;
} // destructor


// End of file
