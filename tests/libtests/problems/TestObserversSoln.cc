// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard == U.S. Geological Survey
// Charles A. Williams == GNS Science
// Matthew G. Knepley == University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California == Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "tests/src/ObserverSolnStub.hh" // USES ObserverSolnStub
#include "tests/src/StubMethodTracker.hh" // USES StubMethodTracker
#include "pylith/problems/ObserversSoln.hh" // USES ObserversSoln
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/problems/problemsfwd.hh" // HOLDSA ObserversSoln

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace problems {
        class TestObserversSoln;
    } // problems
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::problems::TestObserversSoln : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestObserversSoln(void);

    /// Destructor.
    ~TestObserversSoln(void);

    /// Test registerObserver().
    void testRegisterObserver(void);

    /// Test removeObserver().
    void testRemoveObserver(void);

    /// Test setgetTimeScale().
    void testTimeScale(void);

    /// Test verifyObservers().
    void testVerifyObservers(void);

    /// Test notifyObservers().
    void testNotifyObservers(void);

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::problems::ObserversSoln* _observers; ///< Test subject.

    static ObserverSolnStub observerA;
    static ObserverSolnStub observerB;

}; // class TestObserversSoln
pylith::problems::ObserverSolnStub pylith::problems::TestObserversSoln::observerA;
pylith::problems::ObserverSolnStub pylith::problems::TestObserversSoln::observerB;

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestObservesSoln::testRegisterObserver", "[TestObserversSoln]") {
    pylith::problems::TestObserversSoln().testRegisterObserver();
}
TEST_CASE("TestObservesSoln::testRemoveObserver", "[TestObserversSoln]") {
    pylith::problems::TestObserversSoln().testRemoveObserver();
}
TEST_CASE("TestObservesSoln::testTimeScale", "[TestObserversSoln]") {
    pylith::problems::TestObserversSoln().testTimeScale();
}
TEST_CASE("TestObservesSoln::testVerifyObservers", "[TestObserversSoln]") {
    pylith::problems::TestObserversSoln().testVerifyObservers();
}
TEST_CASE("TestObservesSoln::testNotifyObservers", "[TestObserversSoln]") {
    pylith::problems::TestObserversSoln().testNotifyObservers();
}

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::TestObserversSoln::TestObserversSoln(void) :
    _observers(new ObserversSoln()) {
    assert(_observers);
    _observers->registerObserver(&observerA);
    _observers->registerObserver(&observerB);
} // setUp


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::TestObserversSoln::~TestObserversSoln(void) {
    delete _observers;_observers = NULL;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test registerObserver().
void
pylith::problems::TestObserversSoln::testRegisterObserver(void) {
    assert(_observers);
    CHECK(size_t(1) == _observers->_observers.count(&observerA));
    CHECK(size_t(1) == _observers->_observers.count(&observerB));
} // testRegisterObserver


// ------------------------------------------------------------------------------------------------
// Test removeObserver().
void
pylith::problems::TestObserversSoln::testRemoveObserver(void) {
    assert(_observers);
    CHECK(size_t(1) == _observers->_observers.count(&observerA));
    CHECK(size_t(1) == _observers->_observers.count(&observerB));

    _observers->removeObserver(&observerA);
    CHECK(size_t(0) == _observers->_observers.count(&observerA));
    CHECK(size_t(1) == _observers->_observers.count(&observerB));

    _observers->removeObserver(&observerB);
    CHECK(size_t(0) == _observers->_observers.count(&observerA));
    CHECK(size_t(0) == _observers->_observers.count(&observerB));
} // testRegisterObserver


// ------------------------------------------------------------------------------------------------
// Test setgetTimeScale().
void
pylith::problems::TestObserversSoln::testTimeScale(void) {
    assert(_observers);

    // Check default
    PylithReal value = 1.0;
    CHECK(value == observerA.getTimeScale());
    CHECK(value == observerB.getTimeScale());

    // Check set value
    value = 2.0;
    _observers->setTimeScale(value);
    CHECK(value == observerA.getTimeScale());
    CHECK(value == observerB.getTimeScale());
} // testTimeScale


// ------------------------------------------------------------------------------------------------
// Test verifyObservers().
void
pylith::problems::TestObserversSoln::testVerifyObservers(void) {
    assert(_observers);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    pylith::topology::Mesh mesh;
    pylith::topology::Field solution(mesh);
    _observers->verifyObservers(solution);

    CHECK(size_t(2) == tracker.getMethodCount("pylith::problems::ObserverPhysicsStub::verifyConfiguration"));
} // testVerifyObservers


// ------------------------------------------------------------------------------------------------
// Test notifyObservers().
void
pylith::problems::TestObserversSoln::testNotifyObservers(void) {
    assert(_observers);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    const PylithReal t = 1.0;
    const PylithInt tindex = 1;
    pylith::topology::Mesh mesh;
    pylith::topology::Field solution(mesh);
    _observers->notifyObservers(t, tindex, solution);

    CHECK(size_t(2) == tracker.getMethodCount("pylith::problems::ObserverPhysicsStub::update"));
} // testNotifyObservers


// End of file
