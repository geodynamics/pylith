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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "tests/src/ObserverPhysicsStub.hh" // USES ObserverPhysicsStub
#include "tests/src/PhysicsImplementationStub.hh" // USES PhysicsImplementationStub
#include "tests/src/StubMethodTracker.hh" // USES StubMethodTracker
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace problems {
        class TestObserversPhysics;
    } // problems
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::problems::TestObserversPhysics : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestObserversPhysics(void);

    /// Destructor.
    ~TestObserversPhysics(void);

    /// Test registerObserver().
    void testRegisterObserver(void);

    /// Test removeObserver().
    void testRemoveObserver(void);

    /// Test setPhysicsImplemetation().
    void testSetPhysicsImplementation(void);

    /// Test setgetTimeScale().
    void testTimeScale(void);

    /// Test verifyObservers().
    void testVerifyObservers(void);

    /// Test notifyObservers().
    void testNotifyObservers(void);

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::problems::ObserversPhysics* _observers; ///< Test subject.

    static ObserverPhysicsStub observerA;
    static ObserverPhysicsStub observerB;

}; // class TestObserversPhysics

pylith::problems::ObserverPhysicsStub pylith::problems::TestObserversPhysics::observerA;
pylith::problems::ObserverPhysicsStub pylith::problems::TestObserversPhysics::observerB;

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestObservesPhysics::testRegisterObserver", "[TestObserversPhysics]") {
    pylith::problems::TestObserversPhysics().testRegisterObserver();
}
TEST_CASE("TestObservesPhysics::testRemoveObserver", "[TestObserversPhysics]") {
    pylith::problems::TestObserversPhysics().testRemoveObserver();
}
TEST_CASE("TestObservesPhysics::testSetPhysicsImplementation", "[TestObserversPhysics]") {
    pylith::problems::TestObserversPhysics().testSetPhysicsImplementation();
}
TEST_CASE("TestObservesPhysics::testTimeScale", "[TestObserversPhysics]") {
    pylith::problems::TestObserversPhysics().testTimeScale();
}
TEST_CASE("TestObservesPhysics::testVerifyObservers", "[TestObserversPhysics]") {
    pylith::problems::TestObserversPhysics().testVerifyObservers();
}
TEST_CASE("TestObservesPhysics::testNotifyObservers", "[TestObserversPhysics]") {
    pylith::problems::TestObserversPhysics().testNotifyObservers();
}

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::TestObserversPhysics::TestObserversPhysics(void) :
    _observers(new ObserversPhysics()) {
    assert(_observers);
    _observers->registerObserver(&observerA);
    _observers->registerObserver(&observerB);
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::problems::TestObserversPhysics::~TestObserversPhysics(void) {
    delete _observers;_observers = NULL;
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Test registerObserver().
void
pylith::problems::TestObserversPhysics::testRegisterObserver(void) {
    assert(_observers);
    CHECK(size_t(1) == _observers->_observers.count(&observerA));
    CHECK(size_t(1) == _observers->_observers.count(&observerB));
    CHECK(size_t(2) == _observers->size());
} // testRegisterObserver


// ---------------------------------------------------------------------------------------------------------------------
// Test removeObserver().
void
pylith::problems::TestObserversPhysics::testRemoveObserver(void) {
    assert(_observers);
    CHECK(size_t(1) == _observers->_observers.count(&observerA));
    CHECK(size_t(1) == _observers->_observers.count(&observerB));
    CHECK(size_t(2) == _observers->size());

    _observers->removeObserver(&observerA);
    CHECK(size_t(0) == _observers->_observers.count(&observerA));
    CHECK(size_t(1) == _observers->_observers.count(&observerB));
    CHECK(size_t(1) == _observers->size());

    _observers->removeObserver(&observerB);
    CHECK(size_t(0) == _observers->_observers.count(&observerA));
    CHECK(size_t(0) == _observers->_observers.count(&observerB));
    CHECK(size_t(0) == _observers->size());
} // testRegisterObserver


// ---------------------------------------------------------------------------------------------------------------------
// Test setPhysicsImplementation().
void
pylith::problems::TestObserversPhysics::testSetPhysicsImplementation(void) {
    pylith::feassemble::PhysicsImplementationStub physics;

    assert(_observers);
    _observers->setPhysicsImplementation(&physics);
} // testSetPhysicsImplementation


// ---------------------------------------------------------------------------------------------------------------------
// Test setgetTimeScale().
void
pylith::problems::TestObserversPhysics::testTimeScale(void) {
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


// ---------------------------------------------------------------------------------------------------------------------
// Test verifyObservers().
void
pylith::problems::TestObserversPhysics::testVerifyObservers(void) {
    assert(_observers);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    pylith::topology::Mesh mesh;
    pylith::topology::Field solution(mesh);
    _observers->verifyObservers(solution);

    CHECK(size_t(2) == tracker.getMethodCount("pylith::problems::ObserverPhysicsStub::verifyConfiguration"));
} // testVerifyObservers


// ---------------------------------------------------------------------------------------------------------------------
// Test notifyObservers().
void
pylith::problems::TestObserversPhysics::testNotifyObservers(void) {
    assert(_observers);

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    const PylithReal t = 1.0;
    const PylithInt tindex = 1;
    pylith::topology::Mesh mesh;
    pylith::topology::Field solution(mesh);
    _observers->notifyObservers(t, tindex, solution, true);
    _observers->notifyObservers(t, tindex, solution, false);

    CHECK(size_t(4) == tracker.getMethodCount("pylith::problems::ObserverPhysicsStub::update"));
} // testNotifyObservers


// End of file
