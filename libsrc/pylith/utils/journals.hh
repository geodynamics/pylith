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
#define PYLITH_COMPONENT_INFO_ROOT(channel, msg) \
        do { \
            if (pylith::utils::MPI::isRoot()) { \
                pythia::journal::info_t info(channel); \
                info << pythia::journal::at(__HERE__) \
                     << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                     << msg << pythia::journal::endl; } \
        } while (0)

#define PYLITH_COMPONENT_INFO(channel, msg) \
        do { \
            pythia::journal::info_t info(channel); \
            info << pythia::journal::at(__HERE__) \
                 << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                 << msg << pythia::journal::endl; \
        } while (0)

#define PYLITH_COMPONENT_DEBUG(channel, msg) \
        do { \
            pythia::journal::debug_t debug(channel); \
            debug << pythia::journal::at(__HERE__) \
                  << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                  << msg << pythia::journal::endl; \
        } while (0)

#define PYLITH_COMPONENT_WARNING(channel, msg) \
        do { \
            pythia::journal::warning_t warning(channel); \
            warning << pythia::journal::at(__HERE__) \
                    << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                    << msg << pythia::journal::endl; \
        } while (0)

#define PYLITH_COMPONENT_ERROR(ExceptionType, channel, msg) \
        do { \
            pythia::journal::error_t error(channel); \
            error << pythia::journal::at(__HERE__) \
                  << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                  << msg << pythia::journal::endl; \
            throw ExceptionType((pylith::ErrorMessage() << msg)); \
        } while (0)


#define PYLITH_COMPONENT_FIREWALL(ExceptionType, channel, msg) \
        do { \
            pythia::journal::error_t firewall(channel); \
            firewall << pythia::journal::at(__HERE__) \
                     << "Component '"<<PyreComponent::getFullIdentifier()<<"': " \
                     << msg << pythia::journal::endl; \
            throw ExceptionType((pylith::ErrorMessage() << msg)); \
        } while (0)


// General ----------------------------------------------------------------------------------------
#define PYLITH_INFO_ROOT(channel, msg) \
        do { \
            if (pylith::utils::MPI::isRoot()) { \
                pythia::journal::info_t info(channel); \
                info << pythia::journal::at(__HERE__) \
                     << msg << pythia::journal::endl; } \
        } while (0)

#define PYLITH_INFO(channel, msg) \
        do { \
            pythia::journal::info_t info(channel); \
            info << pythia::journal::at(__HERE__) \
                 << msg << pythia::journal::endl; \
        } while (0)

#define PYLITH_DEBUG(channel, msg) \
        do { \
            pythia::journal::debug_t debug(channel); \
            debug << pythia::journal::at(__HERE__) \
                  << msg << pythia::journal::endl; \
        } while (0)

#define PYLITH_WARNING(channel, msg) \
        do { \
            pythia::journal::warning_t warning(channel); \
            warning << pythia::journal::at(__HERE__) \
                    << msg << pythia::journal::endl; \
        } while (0)

#define PYLITH_ERROR(ExceptionType, channel, msg) \
        do { \
            pythia::journal::error_t error(channel); \
            error << pythia::journal::at(__HERE__) \
                  << msg << pythia::journal::endl; \
            throw ExceptionType((pylith::ErrorMessage() << msg)); \
        } while (0)

#define PYLITH_FIREWALL(ExceptionType, channel, msg) \
        do { \
            pythia::journal::error_t firewall(channel); \
            firewall << pythia::journal::at(__HERE__) \
                     << msg << pythia::journal::endl; \
            throw ExceptionType((pylith::ErrorMessage() << msg)); \
        } while (0)


// Channel names ----------------------------------------------------------------------------------
namespace pylith::journal {
    // info
    static const std::string about = "about";
    static const std::string application_flow = "application-flow";
    static const std::string application_flow_detail3 = "application-flow-detail-3";
    static const std::string application_flow_detail5 = "application-flow-detail-5";
    static const std::string debug_parameters = "debug-parameters";

    // warning
    static const std::string user_input = "user-input";
    static const std::string deprecation = "deprecation";

    // error
    // static const std::string user_input = "user-input";
    static const std::string external = "external";
    static const std::string output = "output";

    // firewall}
    static const std::string internal = "internal";
    static const std::string logic = "logic";

    // debug
    // static const std::string application_flow;
    static const std::string auxiliary_fields = "auxiliary-fields";
    static const std::string solution = "solution";
    static const std::string residual = "residual";
    static const std::string mms_test = "mms-test";
    static const std::string full_scale_test = "full-scale-test";
    static const std::string integration_kernels = "integration-kernels";
    static const std::string solver = "solver";
    static const std::string mesh = "mesh";
    static const std::string mesh_full_detail = "mesh detail=5";
} // namespace

// End of file
