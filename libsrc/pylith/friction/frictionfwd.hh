// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file libsrc/friction/frictionfwd.hh
 *
 * @brief Forward declarations for PyLith friction objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_friction_frictionfwd_hh)
#define pylith_friction_frictionfwd_hh

namespace pylith {
    namespace friction {
        class FrictionModel;

        class StaticFriction;
        class SlipWeakening;
        class SlipWeakeningTime;
        class SlipWeakeningTimeStable;
        class RateStateAgeing;
        class TimeWeakening;

    } // friction
} // pylith

#endif // pylith_friction_frictionfwd_hh

// End of file
