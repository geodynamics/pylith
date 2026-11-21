// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


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
    } // namespace petsc
} // namespace pylith


// End of file
