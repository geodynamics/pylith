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

#include "TestPhysics.hh" // Implementation of class methods

#include "pylith/testing/PhysicsStub.hh" // USES PhysicsStub
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory

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
    spatialdata::units::Nondimensional normalizer;
    const PylithReal lengthScale = 3.0;
    normalizer.lengthScale(lengthScale);

    CPPUNIT_ASSERT(_physics);
    _physics->setNormalizer(normalizer);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(lengthScale, _physics->_normalizer->lengthScale(), 1.0e-6);
} // testSetNormalizer


// ---------------------------------------------------------------------------------------------------------------------
// Test registerObserver(), removeObserver(), getObservers().
void
pylith::problems::TestPhysics::testObservers(void) {
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
} // testObservers


// ---------------------------------------------------------------------------------------------------------------------
// Test getKernelConstants().
void
pylith::problems::TestPhysics::testGetKernelConstants(void) {
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad Implement test.", false);
} // testGetKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::problems::TestPhysics::testVerifyConfiguration(void) {
    CPPUNIT_ASSERT(_physics);
    CPPUNIT_ASSERT(_solution);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    _physics->verifyConfiguration(*_solution);
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount("pylith::problems::PhysicsStub::verifyConfiguration"));
} // testVerifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Test createIntegrator().
void
pylith::problems::TestPhysics::testCreateIntegrator(void) {
    CPPUNIT_ASSERT(_physics);
    CPPUNIT_ASSERT(_solution);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    _physics->createIntegrator(*_solution);
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount("pylith::problems::PhysicsStub::createIntegrator"));
} // testCreateIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Test createConstraint().
void
pylith::problems::TestPhysics::testCreateConstraint(void) {
    CPPUNIT_ASSERT(_physics);
    CPPUNIT_ASSERT(_solution);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    _physics->createConstraint(*_solution);
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount("pylith::problems::PhysicsStub::createConstraint"));
} // testCreateConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Test createAuxiliaryField().
void
pylith::problems::TestPhysics::testCreateAuxiliaryField(void) {
    CPPUNIT_ASSERT(_physics);
    CPPUNIT_ASSERT(_solution);
    CPPUNIT_ASSERT(_mesh);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    _physics->createAuxiliaryField(*_solution, *_mesh);
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount("pylith::problems::PhysicsStub::createAuxiliaryField"));
} // testCreateAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Test createDerivedField().
void
pylith::problems::TestPhysics::testCreateDerivedField(void) {
    CPPUNIT_ASSERT(_physics);
    CPPUNIT_ASSERT(_solution);
    CPPUNIT_ASSERT(_mesh);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    _physics->createAuxiliaryField(*_solution, *_mesh);
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount("pylith::problems::PhysicsStub::createAuxiliaryField"));
} // testCreateDerivedField


// End of file
