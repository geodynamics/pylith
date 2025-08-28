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

#include "pylith/utils/types.hh" // HASA PylithScalar

#include <limits>

namespace pylith {
    static const double g_acc = 9.80665;

    static const double max_double = std::numeric_limits<double>::max();
    static const float max_float = std::numeric_limits<float>::max();
    static const PylithInt max_int = PETSC_MAX_INT;
    static const PylithInt min_int = PETSC_MIN_INT;
    static const PylithReal max_real = (sizeof(PylithReal) == sizeof(double)) ? max_double : max_float;

}

// End of file
