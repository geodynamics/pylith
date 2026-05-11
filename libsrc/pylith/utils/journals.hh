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

#include <portinfo>
#include "pythia/journal/diagnostics.h"
#include "pylith/utils/mpi.hh"


// Component --------------------------------------------------------------------------------------
#define PYLITH_COMPONENT_INFO_ROOT_NEW(channel, msg) \
    do { \
        if (pylith::utils::MPI::isRoot()) { \
            pythia::journal::info_t info(channel); \
            info << pythia::journal::at(__HERE__) \
                 << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                 << msg << pythia::journal::endl; } \
    } while (0)

#define PYLITH_COMPONENT_INFO_NEW(channel, msg) \
    do { \
        pythia::journal::info_t info(channel); \
        info << pythia::journal::at(__HERE__) \
             << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
             << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_COMPONENT_DEBUG_NEW(channel, msg) \
    do { \
        pythia::journal::debug_t debug(channel); \
        debug << pythia::journal::at(__HERE__) \
              << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
              << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_COMPONENT_WARNING_NEW(channel, msg) \
    do { \
        pythia::journal::warning_t warning(channel); \
        warning << pythia::journal::at(__HERE__) \
                << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_COMPONENT_ERROR_NEW(ExceptionType, channel, msg) \
    do { \
        pythia::journal::error_t error(channel); \
        error << pythia::journal::at(__HERE__) \
              << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
              << msg << pythia::journal::endl; \
        throw ExceptionType((::pylith::ErrorMessage() << msg)); \
    } while (0)


#define PYLITH_COMPONENT_FIREWALL(ExceptionType, channel, msg) \
    do { \
        pythia::journal::error_t firewall(channel); \
        firewall << pythia::journal::at(__HERE__) \
                 << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                 << msg << pythia::journal::endl; \
        throw ExceptionType((::pylith::ErrorMessage() << msg)); \
    } while (0)


// General ----------------------------------------------------------------------------------------
#define PYLITH_INFO_ROOT_NEW(channel, msg) \
    do { \
        if (pylith::utils::MPI::isRoot()) { \
            pythia::journal::info_t info(channel); \
            info << pythia::journal::at(__HERE__) \
                 << msg << pythia::journal::endl; } \
    } while (0)

#define PYLITH_INFO_NEW(channel, msg) \
    do { \
        pythia::journal::info_t info(channel); \
        info << pythia::journal::at(__HERE__) \
             << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_DEBUG_NEW(channel, msg) \
    do { \
        pythia::journal::debug_t debug(channel); \
        debug << pythia::journal::at(__HERE__) \
              << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_WARNING_NEW(channel, msg) \
    do { \
        pythia::journal::warning_t warning(channel); \
        warning << pythia::journal::at(__HERE__) \
                << msg << pythia::journal::endl; \
    } while (0)

#define PYLITH_ERROR_NEW(ExceptionType, channel, msg) \
    do { \
        pythia::journal::error_t error(channel); \
        error << pythia::journal::at(__HERE__) \
              << msg << pythia::journal::endl; \
        throw ExceptionType((::pylith::ErrorMessage() << msg)); \
    } while (0)

#define PYLITH_FIREWALL(ExceptionType, channel, msg) \
    do { \
        pythia::journal::error_t firewall(channel); \
        firewall << pythia::journal::at(__HERE__) \
                 << msg << pythia::journal::endl; \
        throw ExceptionType((::pylith::ErrorMessage() << msg)); \
    } while (0)


// Obsolete ---------------------------------------------------------------------------------------
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

#define PYLITH_COMPONENT_DEBUG(msg) \
    do { \
        pythia::journal::debug_t debug(PyreComponent::getName()); \
        debug << pythia::journal::at(__HERE__) \
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


#define PYLITH_INTERNAL_ERROR(ExceptionType, channel, msg) \
    do { \
        pythia::journal::error_t firewall(channel); \
        firewall << pythia::journal::at(__HERE__) \
                 << msg << pythia::journal::endl; \
        throw ExceptionType((::pylith::ErrorMessage() << msg)); \
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


// Channel names ----------------------------------------------------------------------------------
namespace pylith::journal {
    // info
    static const std::string about = "about";
    static const std::string application_flow = "application-flow";
    static const std::string debug_config = "debug-config";
    static const std::string initialization = "initialization";

    // warning
    static const std::string user_input = "user-input";
    static const std::string deprecation = "deprecation";

    // error
    static const std::string configuration_error = "configuration";
    static const std::string input_error = "input";
    static const std::string output_error = "output";
    static const std::string external_error = "external";

    // firewall}
    static const std::string internal_error = "internal";
    static const std::string logic_error = "logic";

    // debug
    // static const std::string application_flow;
    static const std::string auxiliary_fields = "auxiliary-fields";
    static const std::string solution_fields = "solution-fields";
    static const std::string integration_kernels = "integration-kernels";
    static const std::string mesh = "mesh";
    static const std::string mms_test = "mms-test";
    static const std::string solver = "solver";
} // namespace

// End of file
