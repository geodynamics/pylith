// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

namespace pylith {
    namespace petsc {
        class Application;
        class ErrorHandler;
        class EventLogger;

        class PetscVersion;
        class Options;
        class Defaults;

        class KSPMonitor;
        class SNESMonitor;
        class TSMonitor;

        class TSAdaptImpulse;
    } // petsc
} // pylith


// End of file
