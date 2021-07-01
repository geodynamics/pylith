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

#include "TestObserversSoln.hh" // Implementation of class methods

#include "pylith/problems/ObserversSoln.hh" // USES ObserversSoln
#include "pylith/testing/ObserverSolnStub.hh" // USES ObserverSolnStub
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/testing/StubMethodTracker.hh" // USES StubMethodTracker

namespace pylith {
    namespace problems {
        class _TestObserversSoln {
public:

            static ObserverSolnStub observerA;
            static ObserverSolnStub observerB;
        }; // _TestObserversSoln
        ObserverSolnStub _TestObserversSoln::observerA;
        ObserverSolnStub _TestObserversSoln::observerB;
    } // problems
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::problems::TestObserversSoln);

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::problems::TestObserversSoln::setUp(void) {
    _observers = new ObserversSoln();CPPUNIT_ASSERT(_observers);
    _observers->registerObserver(&_TestObserversSoln::observerA);
    _observers->registerObserver(&_TestObserversSoln::observerB);
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::problems::TestObserversSoln::tearDown(void) {
    delete _observers;_observers = NULL;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Test registerObserver().
void
pylith::problems::TestObserversSoln::testRegisterObserver(void) {
    CPPUNIT_ASSERT(_observers);
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObserversSoln::observerA));
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObserversSoln::observerB));
} // testRegisterObserver


// ---------------------------------------------------------------------------------------------------------------------
// Test removeObserver().
void
pylith::problems::TestObserversSoln::testRemoveObserver(void) {
    CPPUNIT_ASSERT(_observers);
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObserversSoln::observerA));
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObserversSoln::observerB));

    _observers->removeObserver(&_TestObserversSoln::observerA);
    CPPUNIT_ASSERT_EQUAL(size_t(0), _observers->_observers.count(&_TestObserversSoln::observerA));
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObserversSoln::observerB));

    _observers->removeObserver(&_TestObserversSoln::observerB);
    CPPUNIT_ASSERT_EQUAL(size_t(0), _observers->_observers.count(&_TestObserversSoln::observerA));
    CPPUNIT_ASSERT_EQUAL(size_t(0), _observers->_observers.count(&_TestObserversSoln::observerB));
} // testRegisterObserver


// ---------------------------------------------------------------------------------------------------------------------
// Test setgetTimeScale().
void
pylith::problems::TestObserversSoln::testTimeScale(void) {
    CPPUNIT_ASSERT(_observers);

    // Check default
    PylithReal value = 1.0;
    CPPUNIT_ASSERT_EQUAL(value, _TestObserversSoln::observerA.getTimeScale());
    CPPUNIT_ASSERT_EQUAL(value, _TestObserversSoln::observerB.getTimeScale());

    // Check set value
    value = 2.0;
    _observers->setTimeScale(value);
    CPPUNIT_ASSERT_EQUAL(value, _TestObserversSoln::observerA.getTimeScale());
    CPPUNIT_ASSERT_EQUAL(value, _TestObserversSoln::observerB.getTimeScale());
} // testTimeScale


// ---------------------------------------------------------------------------------------------------------------------
// Test verifyObservers().
void
pylith::problems::TestObserversSoln::testVerifyObservers(void) {
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
pylith::problems::TestObserversSoln::testNotifyObservers(void) {
    CPPUNIT_ASSERT(_observers);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    const PylithReal t = 1.0;
    const PylithInt tindex = 1;
    pylith::topology::Mesh mesh;
    pylith::topology::Field solution(mesh);
    _observers->notifyObservers(t, tindex, solution);

    CPPUNIT_ASSERT_EQUAL(size_t(2), tracker.getMethodCount("pylith::problems::ObserverPhysicsStub::update"));
} // testNotifyObservers


// End of file
