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

#include <string>
#include <sstream>

namespace pylith {
    class ErrorMessage;
}

class pylith::ErrorMessage {
public:

    inline ErrorMessage(const char* text);

    inline ErrorMessage(const std::string& text);

    ErrorMessage() = default;

    /// Append any streamable value.
    template <typename T>
    inline
    ErrorMessage& operator<<(const T& value);

    /// Extract the assembled string.
    inline std::string str() const;

    /// Implicit conversion so ErrorMessage can be passed where std::string
    /// is expected.
    inline operator std::string() const;

private:

    std::ostringstream _buffer;

};

#include "ErrorMessage.icc"
