// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2019 University of California, Davis
//
// See COPYING for license information.
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
    } // problems

} // pylith

#endif // pylith_testing_testingfwd_hh

// End of file
