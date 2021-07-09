// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/** @file libsrc/testing/testingfwd.hh
 *
 * @brief Forward declarations for stubs of PyLith objects used in unit testing.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_testing_testingfwd_hh)
#define pylith_testing_testingfwd_hh

namespace pylith {
    namespace testing {
        class StubMethodTracker;
        class MMSTest;
        class FieldTester;
        class TestDriver;
    } // testing

    namespace faults {
        class FaultCohesiveStub;
    } // faults

    namespace feassemble {
        class PhysicsImplementationStub;
    } // feassemble

    namespace problems {
        class PhysicsStub;
        class ObserverSolnStub;
        class ObserverPhysicsStub;
        class ProgressMonitorStub;
    } // problems

} // pylith

#endif // pylith_testing_testingfwd_hh

// End of file
