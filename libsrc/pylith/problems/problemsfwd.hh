// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file libsrc/problems/problemsfwd.hh
 *
 * @brief Forward declarations for PyLith problems objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_problems_problemsfwd_hh)
#define pylith_problems_problemsfwd_hh

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

#endif // pylith_problems_problemsfwd_hh

// End of file
