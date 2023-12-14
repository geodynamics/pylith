// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

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

// End of file
