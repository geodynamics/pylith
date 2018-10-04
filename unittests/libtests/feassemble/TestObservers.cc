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

#include "TestObservers.hh" // Implementation of class methods

#include "pylith/feassemble/Observers.hh" // USES Observers
#include "pylith/feassemble/ObserverStub.hh" // USES ObserverStub
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

namespace pylith {
    namespace feassemble {
        class _TestObservers {
public:

            static ObserverStub observerA;
            static ObserverStub observerB;
        }; // _TestObservers
        ObserverStub _TestObservers::observerA;
        ObserverStub _TestObservers::observerB;
    } // feassemble
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestObservers::setUp(void) {
    _observers = new Observers();CPPUNIT_ASSERT(_observers);
    _observers->registerObserver(&_TestObservers::observerA);
    _observers->registerObserver(&_TestObservers::observerB);
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::feassemble::TestObservers::tearDown(void) {
    delete _observers;_observers = NULL;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Test registerObserver().
void
pylith::feassemble::TestObservers::testRegisterObserver(void) {
    CPPUNIT_ASSERT(_observers);
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObservers::observerA));
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObservers::observerB));
} // testRegisterObserver


// ---------------------------------------------------------------------------------------------------------------------
// Test removeObserver().
void
pylith::feassemble::TestObservers::testRemoveObserver(void) {
    CPPUNIT_ASSERT(_observers);
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObservers::observerA));
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObservers::observerB));

    _observers->removeObserver(&_TestObservers::observerA);
    CPPUNIT_ASSERT_EQUAL(size_t(0), _observers->_observers.count(&_TestObservers::observerA));
    CPPUNIT_ASSERT_EQUAL(size_t(1), _observers->_observers.count(&_TestObservers::observerB));

    _observers->removeObserver(&_TestObservers::observerB);
    CPPUNIT_ASSERT_EQUAL(size_t(0), _observers->_observers.count(&_TestObservers::observerA));
    CPPUNIT_ASSERT_EQUAL(size_t(0), _observers->_observers.count(&_TestObservers::observerB));
} // testRegisterObserver


// ---------------------------------------------------------------------------------------------------------------------
// Test verifyObservers().
void
pylith::feassemble::TestObservers::testVerifyObservers(void) {
    CPPUNIT_ASSERT(_observers);
    try {
        pylith::topology::Mesh mesh;
        pylith::topology::Field solution(mesh);
        _observers->verifyObservers(solution);
    } catch (ObserverStubException err) {
        CPPUNIT_ASSERT_EQUAL(ObserverStubException::VERIFIED, err.getMethodCalled());
    } // try/catch
} // testVerifyObservers


// ---------------------------------------------------------------------------------------------------------------------
// Test notifyObservers().
void
pylith::feassemble::TestObservers::testNotifyObservers(void) {
    CPPUNIT_ASSERT(_observers);
    try {
        const PylithReal t = 1.0;
        const PylithInt tindex = 1;
        pylith::topology::Mesh mesh;
        pylith::topology::Field solution(mesh);
        _observers->notifyObservers(t, tindex, solution);
    } catch (ObserverStubException err) {
        CPPUNIT_ASSERT_EQUAL(ObserverStubException::UPDATED, err.getMethodCalled());
    } // try/catch
} // testNotifyObservers


// End of file
