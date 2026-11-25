// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


/** @file libsrc/friction/frictionfwd.hh
 *
 * @brief Forward declarations for PyLith friction objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */
#pragma once

namespace pylith {
    namespace friction {
        class FrictionModel;

        class StaticFriction;
        class SlipWeakening;
        class SlipWeakeningTime;
        class SlipWeakeningTimeStable;
        class RateStateAgeing;
        class TimeWeakening;

    } // namespace friction
} // namespace pylith


// End of file
