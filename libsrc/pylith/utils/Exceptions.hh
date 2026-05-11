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

#include "ErrorMessage.hh"

#include <portinfo>
#include <stdexcept>
#include <string>
#include <vector>

namespace pylith {
    class Error;
    class RunTimeError;
    class ValueError;
    class OutOfRangeError;
    class IOError;
    class TopologyError;

    class ExternalError;

    class InternalError;
    class InternalLogicError;
} // namespace

class pylith::Error : public std::runtime_error {
public:

    /** Default constructor.
     *
     * @param[in] message Error message.
     */
    explicit Error(const pylith::ErrorMessage& message);

    ~Error() noexcept override = default;

    /// Full formatted message, including source location when available.
    const char* what() const noexcept override;

protected:

    /// Rebuild _what so that what() reflects the current state.
    void _buildWhat();

    /// Capture the stack frames into _traceback.
    void _captureTraceback();

    /// Formats the traceback as a single multi-line string for printing.
    std::string _formatTraceback() const;

private:

    std::string _what; ///< Cached full message for what()
    std::vector<std::string> _traceback; ///< Stack frames (empty if unavailable)

    static constexpr int MAX_FRAMES = 64;
};


// Error types ====================================================================================

/// General runtime failure (logic is correct; environment caused the failure).
class pylith::RunTimeError : public Error {
public:

    using Error::Error;
};


/// An argument or configuration value is semantically invalid.
class pylith::ValueError : public Error {
public:

    using Error::Error;
};


/// A numeric or index value lies outside its permissible range.
class pylith::OutOfRangeError : public ValueError {
public:

    using ValueError::ValueError;
};


/// File or stream I/O failed (open, read, write, parse, …).
class pylith::IOError : public Error {
public:

    using Error::Error;
};


/// The mesh or its topology is invalid or inconsistent.
class pylith::TopologyError : public Error {
public:

    using Error::Error;
};


/// Error originating in PETSc library.

class pylith::ExternalError : public Error {
public:

    using Error::Error;
};


/// A bug or violated invariant was detected (precondition, assertion, …).
class pylith::InternalError : public Error {
public:

    using Error::Error;
};


/// An internal data structure is corrupt or inconsistent.
class pylith::InternalLogicError : public InternalError {
public:

    using InternalError::InternalError;
};
