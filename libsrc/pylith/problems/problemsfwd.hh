// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

namespace pylith {
    namespace problems {
        class Problem;
        class TimeDependent;
        class GreensFns;

        class SolutionFactory;
        class Observer;
        class ObserversSoln;
        class ObserverSoln;

        class Physics;
        class ObserversPhysics;
        class ObserverPhysics;

        class InitialCondition;
        class InitialConditionDomain;
        class InitialConditionPatch;

        class ProgressMonitor;
        class ProgressMonitorTime;
        class ProgressMonitorStep;

    } // problems
} // pylith

// End of file
