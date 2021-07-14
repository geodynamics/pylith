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

#include "TestPhysics.hh" // Implementation of class methods

#include "pylith/testing/PhysicsStub.hh" // USES PhysicsStub
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/utils/types.hh" // USES PylithReal
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include "pylith/testing/ObserverPhysicsStub.hh" // USES ObserversPhysicsStub
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics

#include "spatialdata/spatialdb/UniformDB.hh" // USES UniformDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/testing/StubMethodTracker.hh" // USES StubMethodTracker

// ---------------------------------------------------------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::problems::TestPhysics);

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::problems::TestPhysics::setUp(void) {
    _physics = new PhysicsStub();CPPUNIT_ASSERT(_physics);

    _mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    _solution = new pylith::topology::Field(*_mesh);CPPUNIT_ASSERT(_solution);
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::problems::TestPhysics::tearDown(void) {
    delete _physics;_physics = NULL;
    delete _solution;_solution = NULL;
    delete _mesh;_mesh = NULL;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Test setNormalizer().
void
pylith::problems::TestPhysics::testSetNormalizer(void) {
    PYLITH_METHOD_BEGIN;

    spatialdata::units::Nondimensional normalizer;
    const PylithReal lengthScale = 3.0;
    normalizer.setLengthScale(lengthScale);

    CPPUNIT_ASSERT(_physics);
    _physics->setNormalizer(normalizer);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(lengthScale, _physics->_normalizer->getLengthScale(), 1.0e-6);

    PYLITH_METHOD_END;
} // testSetNormalizer


// ---------------------------------------------------------------------------------------------------------------------
// Test setAuxiliaryFieldDB().
void
pylith::problems::TestPhysics::testSetAuxiliaryFieldDB(void) {
    PYLITH_METHOD_BEGIN;

    spatialdata::spatialdb::UniformDB db;
    db.setLabel("test db");

    CPPUNIT_ASSERT(_physics);
    _physics->setAuxiliaryFieldDB(&db);

    const pylith::feassemble::AuxiliaryFactory* factory = _physics->_getAuxiliaryFactory();CPPUNIT_ASSERT(factory);
    const spatialdata::spatialdb::SpatialDB* queryDB = factory->getQueryDB();CPPUNIT_ASSERT(queryDB);
    CPPUNIT_ASSERT_EQUAL(db.getLabel(),queryDB->getLabel());

    PYLITH_METHOD_END;
} // testSetAuxiliaryFieldDB


// ---------------------------------------------------------------------------------------------------------------------
// Test setAuxiliarySubfieldDiscretization();
void
pylith::problems::TestPhysics::testSetAuxiliarySubfieldDiscretization(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_physics);
    _physics->setAuxiliarySubfieldDiscretization("fieldA", 1, 2, 3, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true);
    _physics->setAuxiliarySubfieldDiscretization("fieldB", 2, 3, 4, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POINT_SPACE, true);

    const pylith::feassemble::AuxiliaryFactory* factory = _physics->_getAuxiliaryFactory();CPPUNIT_ASSERT(factory);

    const pylith::topology::FieldBase::Discretization& feA = factory->getSubfieldDiscretization("fieldA");
    CPPUNIT_ASSERT_EQUAL(1, feA.basisOrder);
    CPPUNIT_ASSERT_EQUAL(2, feA.quadOrder);
    CPPUNIT_ASSERT_EQUAL(3, feA.dimension);
    CPPUNIT_ASSERT_EQUAL(true, feA.isBasisContinuous);
    CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, feA.feSpace);

    const pylith::topology::FieldBase::Discretization& feB = factory->getSubfieldDiscretization("fieldB");
    CPPUNIT_ASSERT_EQUAL(2, feB.basisOrder);
    CPPUNIT_ASSERT_EQUAL(3, feB.quadOrder);
    CPPUNIT_ASSERT_EQUAL(4, feB.dimension);
    CPPUNIT_ASSERT_EQUAL(false, feB.isBasisContinuous);
    CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POINT_SPACE, feB.feSpace);

    PYLITH_METHOD_END;
} // testSetAuxiliarySubfieldDiscretization


// ---------------------------------------------------------------------------------------------------------------------
// Test registerObserver(), removeObserver(), getObservers().
void
pylith::problems::TestPhysics::testObservers(void) {
    PYLITH_METHOD_BEGIN;

    pylith::problems::ObserverPhysicsStub observerA;
    pylith::problems::ObserverPhysicsStub observerB;

    CPPUNIT_ASSERT(_physics);
    CPPUNIT_ASSERT_EQUAL(size_t(0), _physics->getObservers()->size());

    _physics->registerObserver(&observerA);
    CPPUNIT_ASSERT_EQUAL(size_t(1), _physics->getObservers()->size());

    _physics->registerObserver(&observerB);
    CPPUNIT_ASSERT_EQUAL(size_t(2), _physics->getObservers()->size());

    _physics->removeObserver(&observerB);
    CPPUNIT_ASSERT_EQUAL(size_t(1), _physics->getObservers()->size());

    _physics->removeObserver(&observerA);
    CPPUNIT_ASSERT_EQUAL(size_t(0), _physics->getObservers()->size());

    _physics->removeObserver(&observerA);
    CPPUNIT_ASSERT_EQUAL(size_t(0), _physics->getObservers()->size());

    PYLITH_METHOD_END;
} // testObservers


// ---------------------------------------------------------------------------------------------------------------------
// Test getKernelConstants().
void
pylith::problems::TestPhysics::testGetKernelConstants(void) {
    PYLITH_METHOD_BEGIN;

    const size_t numConstants = 2;
    const PylithReal constantsE[numConstants] = { -1.1, 4.4 };

    CPPUNIT_ASSERT(_physics);
    _physics->_kernelConstants = pylith::real_array(constantsE, numConstants);

    const PylithReal dt = 2.0;
    const pylith::real_array& constants = _physics->getKernelConstants(dt);

    CPPUNIT_ASSERT_EQUAL(numConstants, constants.size());
    for (size_t i = 0; i < numConstants; ++i) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(constantsE[i], constants[i], 1.0e-6);
    } // for

    PYLITH_METHOD_END;
} // testGetKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::problems::TestPhysics::testVerifyConfiguration(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_physics);
    CPPUNIT_ASSERT(_solution);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    _physics->verifyConfiguration(*_solution);
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount("pylith::problems::PhysicsStub::verifyConfiguration"));

    PYLITH_METHOD_END;
} // testVerifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Test createIntegrator().
void
pylith::problems::TestPhysics::testCreateIntegrator(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_physics);
    CPPUNIT_ASSERT(_solution);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    _physics->createIntegrator(*_solution);
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount("pylith::problems::PhysicsStub::createIntegrator"));

    PYLITH_METHOD_END;
} // testCreateIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Test createConstraint().
void
pylith::problems::TestPhysics::testCreateConstraint(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_physics);
    CPPUNIT_ASSERT(_solution);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    _physics->createConstraint(*_solution);
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount("pylith::problems::PhysicsStub::createConstraint"));

    PYLITH_METHOD_END;
} // testCreateConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Test createAuxiliaryField().
void
pylith::problems::TestPhysics::testCreateAuxiliaryField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_physics);
    CPPUNIT_ASSERT(_solution);
    CPPUNIT_ASSERT(_mesh);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    _physics->createAuxiliaryField(*_solution, *_mesh);
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount("pylith::problems::PhysicsStub::createAuxiliaryField"));

    PYLITH_METHOD_END;
} // testCreateAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Test createDerivedField().
void
pylith::problems::TestPhysics::testCreateDerivedField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_physics);
    CPPUNIT_ASSERT(_solution);
    CPPUNIT_ASSERT(_mesh);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    _physics->createAuxiliaryField(*_solution, *_mesh);
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount("pylith::problems::PhysicsStub::createAuxiliaryField"));

    PYLITH_METHOD_END;
} // testCreateDerivedField


// End of file
