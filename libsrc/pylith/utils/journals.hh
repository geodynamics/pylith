// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
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

#include "pythia/journal/diagnostics.h"

#define PYLITH_COMPONENT_DEBUG(msg) \
    do { \
        pythia::journal::debug_t debug(PyreComponent::getName()); \
        debug << pythia::journal::at(__HERE__) \
              << "Component '"<<PyreComponent::getIdentifier()<<"': " \
              << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_COMPONENT_INFO(msg) \
    do { \
        pythia::journal::info_t info(PyreComponent::getName()); \
        info << pythia::journal::at(__HERE__) \
             << "Component '"<<PyreComponent::getIdentifier()<<"': " \
             << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_COMPONENT_WARNING(msg) \
    do { \
        pythia::journal::warning_t warning(PyreComponent::getName()); \
        warning << pythia::journal::at(__HERE__) \
                << "Component '"<<PyreComponent::getIdentifier()<<"': " \
                << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_COMPONENT_ERROR(msg) \
    do { \
        pythia::journal::error_t error(PyreComponent::getName()); \
        error << pythia::journal::at(__HERE__) \
              << "Component '"<<PyreComponent::getIdentifier()<<"': " \
              << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_COMPONENT_LOGICERROR(msg) \
    do { \
        std::ostringstream firewall; \
        firewall << pythia::journal::at(__HERE__) \
                 << "Component '"<<PyreComponent::getIdentifier()<<"': " \
                 << msg; \
        throw std::logic_error(firewall.str().c_str()); \
    } while (0)

#define PYLITH_JOURNAL_DEBUG(msg) \
    do { \
        pythia::journal::debug_t debug(GenericComponent::getName()); \
        debug << pythia::journal::at(__HERE__) \
              << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_JOURNAL_INFO(msg) \
    do { \
        pythia::journal::info_t info(GenericComponent::getName()); \
        info << pythia::journal::at(__HERE__) \
             << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_JOURNAL_WARNING(msg) \
    do { \
        pythia::journal::warning_t warning(GenericComponent::getName()); \
        warning << pythia::journal::at(__HERE__) \
                << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_JOURNAL_ERROR(msg) \
    do { \
        pythia::journal::error_t error(GenericComponent::getName()); \
        error << pythia::journal::at(__HERE__) \
              << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_JOURNAL_LOGICERROR(msg) \
    do { \
        std::ostringstream firewall; \
        firewall << pythia::journal::at(__HERE__) \
                 << msg; \
        throw std::logic_error(firewall.str().c_str()); \
    } while (0)

#endif // pylith_utils_journals_hh

// End of file
