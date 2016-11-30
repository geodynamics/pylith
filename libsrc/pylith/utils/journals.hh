// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/utils/journals.hh
 *
 * @brief Macros for easy use of JournalingComponent.
 */

#if !defined(pylith_utils_journals_hh)
#define pylith_utils_journals_hh

#include "journal/debug.h"
#include "journal/info.h"
#include "journal/warning.h"
#include "journal/error.h"

#include <cassert> // USES assert()


#define PYLITH_JOURNAL_DEBUG(msg) \
    do { \
        if (!_debug) { \
            JournalingComponent::initialize(); assert(_debug); \
        } \
        *_debug << journal::at(__HERE__) << msg << journal::endl; \
    } while(0)

#define PYLITH_JOURNAL_INFO(msg) \
    do { \
        if (!_info) { \
            JournalingComponent::initialize(); assert(_info); \
        } \
        *_info << journal::at(__HERE__) << msg << journal::endl; \
    } while(0)

#define PYLITH_JOURNAL_WARNING(msg) \
    do { \
        if (!_warning) { \
            JournalingComponent::initialize(); assert(_warning); \
        } \
        *_warning << journal::at(__HERE__) << msg << journal::endl; \
    } while(0)

#define PYLITH_JOURNAL_ERROR(msg) \
    do { \
        if (!_error) { \
            JournalingComponent::initialize(); assert(_error); \
        } \
        *_error << journal::at(__HERE__) << msg << journal::endl; \
    } while(0)

#endif // pylith_utils_journals_hh

// End of file
