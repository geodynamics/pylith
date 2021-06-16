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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

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
