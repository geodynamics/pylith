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

#include "TestObserversPhysics.hh" // Implementation of class methods

#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/testing/ObserverPhysicsStub.hh" // USES ObserverPhysicsStub
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/testing/PhysicsImplementationStub.hh" // USES PhysicsImplementationStub
#include "pylith/testing/StubMethodTracker.hh" // USES StubMethodTracker

namespace pylith {
    namespace problems {
        class _TestObserversPhysics {
public:

            static ObserverPhysicsStub observerA;
            static ObserverPhysicsStub observerB;
        }; // _TestObserversPhysics
        ObserverPhysicsStub _TestObserversPhysics::observerA;
        ObserverPhysicsStub _TestObserversPhysics::observerB;
    } // problems
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::problems::TestObserversPhysics);

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::problems::TestObserversPhysics::setUp(void) {
    _observers = new ObserversPhysics();CPPUNIT_ASSERT(_observers);
    _observers->registerObserver(&_TestObserversPhysics::observerA);
    _observers->registerObserver(&_TestObserversPhysics::observerB);
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::problems::TestObserversPhysics::tearDown(void) {
    delete _observers;_observers = NULL;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Test registerObserver().
void
pylith::problems::TestObserversPhysics::testRegisterObserver(void) {
    CPPUNIT_ASSERT(_observers);
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObserversPhysics::observerA));
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObserversPhysics::observerB));
    CPPUNIT_ASSERT_EQUAL(size_t(2), _observers->size());
} // testRegisterObserver


// ---------------------------------------------------------------------------------------------------------------------
// Test removeObserver().
void
pylith::problems::TestObserversPhysics::testRemoveObserver(void) {
    CPPUNIT_ASSERT(_observers);
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObserversPhysics::observerA));
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObserversPhysics::observerB));
    CPPUNIT_ASSERT_EQUAL(size_t(2), _observers->size());

    _observers->removeObserver(&_TestObserversPhysics::observerA);
    CPPUNIT_ASSERT_EQUAL(size_t(0), _observers->_observers.count(&_TestObserversPhysics::observerA));
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObserversPhysics::observerB));
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->size());

    _observers->removeObserver(&_TestObserversPhysics::observerB);
    CPPUNIT_ASSERT_EQUAL(size_t(0), _observers->_observers.count(&_TestObserversPhysics::observerA));
    CPPUNIT_ASSERT_EQUAL(size_t(0), _observers->_observers.count(&_TestObserversPhysics::observerB));
    CPPUNIT_ASSERT_EQUAL(size_t(0), _observers->size());
} // testRegisterObserver


// ---------------------------------------------------------------------------------------------------------------------
// Test setPhysicsImplementation().
void
pylith::problems::TestObserversPhysics::testSetPhysicsImplementation(void) {
    pylith::feassemble::PhysicsImplementationStub physics;

    CPPUNIT_ASSERT(_observers);
    _observers->setPhysicsImplementation(&physics);
} // testSetPhysicsImplementation


// ---------------------------------------------------------------------------------------------------------------------
// Test setgetTimeScale().
void
pylith::problems::TestObserversPhysics::testTimeScale(void) {
    CPPUNIT_ASSERT(_observers);

    // Check default
    PylithReal value = 1.0;
    CPPUNIT_ASSERT_EQUAL(value, _TestObserversPhysics::observerA.getTimeScale());
    CPPUNIT_ASSERT_EQUAL(value, _TestObserversPhysics::observerB.getTimeScale());

    // Check set value
    value = 2.0;
    _observers->setTimeScale(value);
    CPPUNIT_ASSERT_EQUAL(value, _TestObserversPhysics::observerA.getTimeScale());
    CPPUNIT_ASSERT_EQUAL(value, _TestObserversPhysics::observerB.getTimeScale());
} // testTimeScale


// ---------------------------------------------------------------------------------------------------------------------
// Test verifyObservers().
void
pylith::problems::TestObserversPhysics::testVerifyObservers(void) {
    CPPUNIT_ASSERT(_observers);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    pylith::topology::Mesh mesh;
    pylith::topology::Field solution(mesh);
    _observers->verifyObservers(solution);

    CPPUNIT_ASSERT_EQUAL(size_t(2), tracker.getMethodCount("pylith::problems::ObserverPhysicsStub::verifyConfiguration"));
} // testVerifyObservers


// ---------------------------------------------------------------------------------------------------------------------
// Test notifyObservers().
void
pylith::problems::TestObserversPhysics::testNotifyObservers(void) {
    CPPUNIT_ASSERT(_observers);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    const PylithReal t = 1.0;
    const PylithInt tindex = 1;
    pylith::topology::Mesh mesh;
    pylith::topology::Field solution(mesh);
    _observers->notifyObservers(t, tindex, solution, true);
    _observers->notifyObservers(t, tindex, solution, false);

    CPPUNIT_ASSERT_EQUAL(size_t(4), tracker.getMethodCount("pylith::problems::ObserverPhysicsStub::update"));
} // testNotifyObservers


// End of file
