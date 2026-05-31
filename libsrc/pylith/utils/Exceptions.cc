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

#if defined(HAVE_BACKTRACE)
#include <execinfo.h>
#include <cstdlib> // USES std::free
#endif

#include <sstream> // USES std::ostringstream
#include <cstring> // USES std::strlen

// ------------------------------------------------------------------------------------------------
pylith::Error::Error(const ErrorMessage& message)
    : std::runtime_error(message.str()) {
    _captureTraceback();
    _buildWhat();
} // constructor


// ------------------------------------------------------------------------------------------------
const char*
pylith::Error::what() const noexcept {
    return _what.c_str();
} // destructor


// ------------------------------------------------------------------------------------------------
// Add context to existing error message.
void
pylith::Error::addContext(const pylith::ErrorMessage& context) {
    _buildWhat(context.str());
}


// ------------------------------------------------------------------------------------------------
void
pylith::Error::_buildWhat(const std::string& context) {
    std::ostringstream oss;
    oss << std::runtime_error::what();
    if (context.length() > 0) {
        oss << "\n" << context;
    } // if
    oss << _formatTraceback();

    _what = oss.str();
} // _buildWhat


// ------------------------------------------------------------------------------------------------
void
pylith::Error::_captureTraceback() {
#if defined(HAVE_BACKTRACE)
    void*  frames[MAX_FRAMES];
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
#endif
} // captureTraceback


// ------------------------------------------------------------------------------------------------
std::string
pylith::Error::_formatTraceback() const {
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
