// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include "Exceptions.hh"

#if defined(HAVE_CPPTRACE)
#include "cpptrace/cpptrace.hpp"
#else
#if defined(HAVE_BACKTRACE)
#include <execinfo.h>
#include <cstdlib> // USES std::free
#endif
#endif

#include <sstream> // USES std::ostringstream
#include <cstring> // USES std::strlen

// ------------------------------------------------------------------------------------------------
pylith::exceptions::Error::Error(const pylith::exceptions::ErrorMessage& message,
                                 const char* filename,
                                 const long line,
                                 const char* function) :
    std::runtime_error(message.str()),
    _line(line),
    _function(function),
    _filename(filename) {
    _captureTraceback();
    _buildWhat();
} // constructor


// ------------------------------------------------------------------------------------------------
const char*
pylith::exceptions::Error::what() const noexcept {
    return _what.c_str();
} // destructor


// ------------------------------------------------------------------------------------------------
// Add context to existing error message.
void
pylith::exceptions::Error::addContext(const pylith::exceptions::ErrorMessage& context) {
    _buildWhat(context.str());
}


// ------------------------------------------------------------------------------------------------
void
pylith::exceptions::Error::_buildWhat(const std::string& context) {
    std::ostringstream oss;
    oss << std::runtime_error::what() << "\n"
        << "    File: \"" << _filename << "\", line " << _line << "\n"
        << "    Function: " << _function << "\n";
    if (context.length() > 0) {
        oss << "\n" << context;
    } // if
    oss << _formatTraceback();

    _what = oss.str();
} // _buildWhat


// ------------------------------------------------------------------------------------------------
void
pylith::exceptions::Error::_captureTraceback() {
    size_t numFrames = 0;
#if defined(HAVE_CPPTRACE)
    const size_t skipFrames = 6; // Skip exception lines in traceback

    cpptrace::stacktrace trace = cpptrace::generate_trace(skipFrames);

    for (const auto& frame : trace.frames) {
        std::ostringstream oss;
        oss << (frame.symbol.empty() ? "??" : frame.symbol)
            << " at "
            << (frame.filename.empty() ? "??" : frame.filename)
            << ":"
            << frame.line.value_or(0);

        _traceback.push_back(oss.str());
    } // for
    numFrames = trace.frames.size();

#else
#if defined(HAVE_BACKTRACE)
    void* frames[MAX_FRAMES];
    const size_t skipFrames = 3; // Skip exception lines in traceback

    size_t numFrames = ::backtrace(frames, MAX_FRAMES);
    if (numFrames <= 0) {return;}

    char** symbols = ::backtrace_symbols(frames, numFrames);
    if (!symbols) {return;}

    _traceback.reserve(numFrames);
    for (size_t i = skipFrames; i < numFrames; ++i) {
        _traceback.emplace_back(symbols[i] ? symbols[i] : "<unknown>");
    }
    std::free(symbols);
    numFrames -= skipFrames;

#endif
#endif

    // Remove libpython from traceback
    size_t last = 0;
    for (size_t i = 0; i < numFrames; ++i) {
        const size_t index = numFrames - 1 - i;
        if (_traceback[index].find("libpython") != std::string::npos) {
            last = index;
            break;
        } // if
    } // for
    if (last > 0) {
        size_t first = last;
        for (size_t i = 0; i < last; ++i) {
            const size_t index = last - 1 - i;
            if (_traceback[index].find("libpython") != std::string::npos) {
                first = index;
            } else {
                break;
            } // if/else
        } // for
        if (last > first) {
            size_t newSize = numFrames - (last-first) - 1;
            for (size_t i = first; i < newSize; ++i) {
                _traceback[i] = _traceback[last + 1 + (i-first)];
            } // for
            _traceback.resize(newSize);
        } // if
    } // if
} // captureTraceback


// ------------------------------------------------------------------------------------------------
std::string
pylith::exceptions::Error::_formatTraceback() const {
    if (_traceback.empty()) {
        return "(backtrace not available)\n";
    } // if

    std::ostringstream oss;
    oss << "\nC++ traceback (" << _traceback.size() << " frames):\n";
    for (std::size_t i = 0; i < _traceback.size(); ++i) {
        oss << "  [" << i << "]  " << _traceback[i] << "\n";
    }
    return oss.str();
} // formatTraceback
