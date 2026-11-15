// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 All rights reserved.
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
    } // namespace testing

    namespace faults {
        class FaultCohesiveStub;
    } // namespace faults

    namespace feassemble {
        class PhysicsImplementationStub;
    } // namespace feassemble

    namespace problems {
        class PhysicsStub;
        class ObserverSolnStub;
        class ObserverPhysicsStub;
        class ProgressMonitorStub;
    } // namespace problems

} // namespace pylith

// End of file
