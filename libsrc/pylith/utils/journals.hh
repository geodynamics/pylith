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

#include "pythia/journal/diagnostics.h"
#include "pylith/utils/mpi.hh"

#define PYLITH_COMPONENT_DEBUG(msg) \
        do { \
            pythia::journal::debug_t debug(PyreComponent::getName()); \
            debug << pythia::journal::at(__HERE__) \
                  << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                  << msg << pythia::journal::endl; \
        } while (0)

#define PYLITH_COMPONENT_INFO_ROOT(msg) \
        do { \
            if (pylith::utils::MPI::isRoot()) { \
                pythia::journal::info_t info(PyreComponent::getName()); \
                info << pythia::journal::at(__HERE__) \
                     << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                     << msg << pythia::journal::endl; } \
        } while (0)

#define PYLITH_COMPONENT_INFO(msg) \
        do { \
            pythia::journal::info_t info(PyreComponent::getName()); \
            info << pythia::journal::at(__HERE__) \
                 << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                 << msg << pythia::journal::endl; \
        } while (0)

#define PYLITH_COMPONENT_WARNING(msg) \
        do { \
            pythia::journal::warning_t warning(PyreComponent::getName()); \
            warning << pythia::journal::at(__HERE__) \
                    << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                    << msg << pythia::journal::endl; \
        } while (0)

#define PYLITH_COMPONENT_ERROR(msg) \
        do { \
            pythia::journal::error_t error(PyreComponent::getName()); \
            error << pythia::journal::at(__HERE__) \
                  << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                  << msg << pythia::journal::endl; \
        } while (0)

#define PYLITH_COMPONENT_LOGICERROR(msg) \
        do { \
            std::ostringstream firewall; \
            firewall << pythia::journal::at(__HERE__) \
                     << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                     << msg; \
            throw std::logic_error(firewall.str().c_str()); \
        } while (0)

#define PYLITH_JOURNAL_DEBUG(msg) \
        do { \
            pythia::journal::debug_t debug(GenericComponent::getName()); \
            debug << pythia::journal::at(__HERE__) \
                  << msg << pythia::journal::endl; \
        } while (0)

#define PYLITH_JOURNAL_INFO_ROOT(msg) \
        do { \
            if (pylith::utils::MPI::isRoot()) { \
                pythia::journal::info_t info(GenericComponent::getName()); \
                info << pythia::journal::at(__HERE__) \
                     << msg << pythia::journal::endl; } \
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

// End of file
