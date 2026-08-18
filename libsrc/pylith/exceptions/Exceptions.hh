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
    namespace exceptions {
        class Error;
        class RuntimeError;

        class ValueError;
        class InvalidParameterError;
        class EmptyStringError;
        class InvalidUnitVectorError;
        class SubfieldNotFoundError;
        class LabelNotFoundError;
        class InvalidNameError;
        class OutOfRangeError;

        class IOError;

        class TopologyError;

        class ExternalError;
        class PetscError;

        class InternalError;
        class InternalLogicError;
        class SwitchLogicError;
        class NotImplementedError;
    } // namespace
} // namespace

class pylith::exceptions::Error : public std::runtime_error {
public:

    /** Default constructor.
     *
     * @param[in] message Error message.
     */
    explicit Error(const pylith::exceptions::ErrorMessage& message,
                   const char* filename,
                   const long line,
                   const char* function);

    ~Error() noexcept override = default;

    /// Full formatted message, including source location when available.
    const char* what() const noexcept override;

    /// Add context to existing error message.
    void addContext(const pylith::exceptions::ErrorMessage& context);

protected:

    /// Rebuild _what so that what() reflects the current state.
    void _buildWhat(const std::string& context="");

    /// Capture the stack frames into _traceback.
    void _captureTraceback();

    /// Formats the traceback as a single multi-line string for printing.
    std::string _formatTraceback() const;

private:

    int _line;
    std::string _function;
    std::string _filename;

    std::string _what; ///< Cached full message for what()
    std::vector<std::string> _traceback; ///< Stack frames (empty if unavailable)

    static constexpr int MAX_FRAMES = 64;
};


// Error types ====================================================================================

/// General runtime failure (logic is correct; environment caused the failure).
class pylith::exceptions::RuntimeError : public Error {
public:

    using Error::Error;
};


// Value errors ===================================================================================

/// An argument or configuration value is semantically invalid.
class pylith::exceptions::ValueError : public Error {
public:

    using Error::Error;
};

/// Invalid parameter choice.
class pylith::exceptions::InvalidParameterError : public ValueError {
public:

    using ValueError::ValueError;
};

/// Empty string is not valid.
class pylith::exceptions::EmptyStringError : public ValueError {
public:

    using ValueError::ValueError;
};

/// Invalid unit vector.
class pylith::exceptions::InvalidUnitVectorError : public ValueError {
public:

    using ValueError::ValueError;
};

/// Subfield not found.
class pylith::exceptions::SubfieldNotFoundError : public ValueError {
public:

    using ValueError::ValueError;
};

/// Label not found.
class pylith::exceptions::LabelNotFoundError : public ValueError {
public:

    using ValueError::ValueError;
};

/// Label not found.
class pylith::exceptions::InvalidNameError : public ValueError {
public:

    using ValueError::ValueError;
};

/// A numeric or index value lies outside its permissible range.
class pylith::exceptions::OutOfRangeError : public ValueError {
public:

    using ValueError::ValueError;
};


/// File or stream I/O failed (open, read, write, parse, …).
class pylith::exceptions::IOError : public Error {
public:

    using Error::Error;
};


/// The mesh or its topology is invalid or inconsistent.
class pylith::exceptions::TopologyError : public Error {
public:

    using Error::Error;
};


/// Error originating in PETSc library.

class pylith::exceptions::ExternalError : public Error {
public:

    using Error::Error;
};


class pylith::exceptions::PetscError : public Error {
public:

    using Error::Error;
};


/// A bug or violated invariant was detected (precondition, assertion, …).
class pylith::exceptions::InternalError : public Error {
public:

    using Error::Error;
};


/// An internal data structure is corrupt or inconsistent.
class pylith::exceptions::InternalLogicError : public InternalError {
public:

    using InternalError::InternalError;
};

/// Switch statement logic error (unknown case).
class pylith::exceptions::SwitchLogicError : public InternalLogicError {
public:

    using InternalLogicError::InternalLogicError;
};

/// Not implemented.
class pylith::exceptions::NotImplementedError : public InternalLogicError {
public:

    using InternalLogicError::InternalLogicError;
};
