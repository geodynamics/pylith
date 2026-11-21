// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


#pragma once

#include <sstream>
#include <string>

namespace pylith {
    namespace exceptions {
        class ErrorMessage;
    }
} // namespace pylith

class pylith::exceptions::ErrorMessage {
  public:
    inline ErrorMessage(const char* text);

    inline ErrorMessage(const std::string& text);

    ErrorMessage() = default;

    /// Append any streamable value.
    template <typename T> inline ErrorMessage& operator<<(const T& value);

    /// Extract the assembled string.
    inline std::string str() const;

    /// Implicit conversion so ErrorMessage can be passed where std::string
    /// is expected.
    inline operator std::string() const;

  private:
    std::ostringstream _buffer;
};

#include "ErrorMessage.icc"


// End of file
