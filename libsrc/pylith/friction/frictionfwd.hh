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
