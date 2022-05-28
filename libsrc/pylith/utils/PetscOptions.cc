// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "PetscOptions.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/materials/Material.hh" // USES Material

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

namespace pylith {
    namespace utils {
        class _PetscOptions {
public:

            /// Write options to journal;
            static
            void write(pythia::journal::info_t& info,
                       const char* heading,
                       const PetscOptions& options);

            /** Check if simulation is running in parallel.
             *
             * @param[in] solution Solution field for problem.
             * @returns True if solving problem in parallel, False if in serial.
             */
            static
            bool isParallel(const pylith::topology::Field& solution);

            /** Check if simulation is has a fault.
             *
             * @param[in] solution Solution field for problem.
             * @returns True if problem has a fault, false otherwise;
             */
            static
            bool hasFault(const pylith::topology::Field& solution);

            /** Add debugging options.
             *
             * @param[in] options PETSc options.
             */
            static
            void addTesting(PetscOptions* options);

            /** Add monitoring options.
             *
             * @param[in] options PETSc options.
             */
            static
            void addMonitoring(PetscOptions* options);

            /** Add default solver tolerances to options.
             *
             * @param[in] options PETSc options.
             */
            static
            void addSolverTolerances(PetscOptions* options);

        };
    }
}

// ------------------------------------------------------------------------------------------------
const int pylith::utils::PetscDefaults::NONE = 0x0;
const int pylith::utils::PetscDefaults::MONITORS = 0x1;
const int pylith::utils::PetscDefaults::SOLVER = 0x2;
const int pylith::utils::PetscDefaults::PARALLEL = 0x4;
const int pylith::utils::PetscDefaults::TESTING = 0x8;

// ------------------------------------------------------------------------------------------------
// Set default PETSc solver options based on solution field and material.
void
pylith::utils::PetscDefaults::set(const pylith::topology::Field& solution,
                                  const pylith::materials::Material* material,
                                  const int flags) {
    PYLITH_METHOD_BEGIN;
    assert(material);

    if (!flags) {
        PYLITH_METHOD_END;
    } // if

    PetscOptions* options = NULL;
    if (flags & SOLVER) {
        const bool isParallel = flags & PARALLEL || _PetscOptions::isParallel(solution);
        const bool hasFault = _PetscOptions::hasFault(solution);
        options = material->getSolverDefaults(isParallel, hasFault);
    } // if
    if (!options) {
        options = new PetscOptions();
    } // if
    assert(options);

    _PetscOptions::addSolverTolerances(options);
    if (flags & TESTING) {
        _PetscOptions::addTesting(options);
    } // if
    if (flags & MONITORS) {
        _PetscOptions::addMonitoring(options);
    } // if

    options->set();
    delete options;options = NULL;

    PYLITH_METHOD_END;
} // setDefaults


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::utils::PetscOptions::PetscOptions(void) {
    GenericComponent::setName("petscoptions");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::utils::PetscOptions::~PetscOptions(void) {}


// ------------------------------------------------------------------------------------------------
// Add PETSc option.
void
pylith::utils::PetscOptions::add(const char* name,
                                 const char* value) {
    _options[std::string(name)] = value;
} // add


// ------------------------------------------------------------------------------------------------
// Remove PETSc option.
void
pylith::utils::PetscOptions::remove(const char* name) {
    const options_t::iterator iter = _options.find(std::string(name));
    if (_options.end() != iter) {
        _options.erase(iter);
    } // if
} // remove


// ------------------------------------------------------------------------------------------------
// Clear PETSc options.
void
pylith::utils::PetscOptions::clear(void) {
    _options.clear();
} // clear


// ------------------------------------------------------------------------------------------------
// Set PETSc options.
void
pylith::utils::PetscOptions::set(void) {
    PYLITH_METHOD_BEGIN;
    PetscErrorCode err = 0;

    PetscOptions optionsUsed;
    PetscOptions optionsIgnored;
    for (options_t::iterator iter = _options.begin(); iter != _options.end(); ++iter) {
        const char* name = iter->first.c_str();
        const char* value = iter->second.c_str();

        PetscBool exists = PETSC_FALSE;
        err = PetscOptionsHasName(NULL, NULL, name, &exists);PYLITH_CHECK_ERROR(err);
        if (!exists) {
            err = PetscOptionsSetValue(NULL, name, value);PYLITH_CHECK_ERROR(err);
            optionsUsed.add(name, value);
        } else {
            optionsIgnored.add(name, value);
        } // if/else
    } // for

    pythia::journal::info_t info(GenericComponent::getName());
    if (info.state()) {
        _PetscOptions::write(info, "Setting PETSc options:", optionsUsed);
        if (optionsIgnored._options.size() > 0) {
            _PetscOptions::write(info, "Ignoring PETSc options (already set):", optionsIgnored);
        } // if
    } // if

    PYLITH_METHOD_END;
} // set


// ------------------------------------------------------------------------------------------------
// Set PETSc options, overriding any previously set options with the same name.
void
pylith::utils::PetscOptions::override (void) {
    PYLITH_METHOD_BEGIN;
    PetscErrorCode err = 0;

    for (options_t::iterator iter = _options.begin(); iter != _options.end(); ++iter) {
        const char* name = iter->first.c_str();
        const char* value = iter->second.c_str();

        err = PetscOptionsSetValue(NULL, name, value);PYLITH_CHECK_ERROR(err);
    } // for

    pythia::journal::info_t info(GenericComponent::getName());
    if (info.state()) {
        _PetscOptions::write(info, "Setting PETSc options:", (*this));
    } // if

    PYLITH_METHOD_END;
} // set

// ------------------------------------------------------------------------------------------------
// Write options to journal;
void
pylith::utils::_PetscOptions::write(pythia::journal::info_t& info,
                                    const char* heading,
                                    const PetscOptions& options) {
    PYLITH_METHOD_BEGIN;

    info << pythia::journal::at(__HERE__)
         << heading << "\n";
    const PetscOptions::options_t::const_iterator begin = options._options.begin();
    const PetscOptions::options_t::const_iterator end = options._options.end();
    for (PetscOptions::options_t::const_iterator iter = begin; iter != end; ++iter) {
        const std::string name = iter->first.substr(1);
        const std::string value = iter->second;
        if (iter->second.empty()) {
            info << name << " = true\n";
        } else {
            info << name << " = " << value << "\n";
        } // if/else
    } // for
    info << pythia::journal::endl;

    PYLITH_METHOD_END;
} // write


// ------------------------------------------------------------------------------------------------
// Check if simulation is running in parallel.
bool
pylith::utils::_PetscOptions::isParallel(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;

    MPI_Comm comm = solution.getMesh().getComm();
    int numProcs = 0;
    MPI_Comm_size(comm, &numProcs);

    PYLITH_METHOD_RETURN(numProcs > 1);
} // isParallel


// ------------------------------------------------------------------------------------------------
// Check if simulation is has a fault.
bool
pylith::utils::_PetscOptions::hasFault(const pylith::topology::Field& solution) {
    return solution.hasSubfield("lagrange_multiplier_fault");
} // hasFault


// ------------------------------------------------------------------------------------------------
// Add debugging options.
void
pylith::utils::_PetscOptions::addTesting(PetscOptions* options) {
    assert(options);

    // -checkstack only works with PetscInitialize()
    options->add("-malloc_dump");
} // setDebugging


// ------------------------------------------------------------------------------------------------
// Add monitoring options.
void
pylith::utils::_PetscOptions::addMonitoring(PetscOptions* options) {
    assert(options);

    options->add("-ksp_converged_reason");

    options->add("-snes_converged_reason");
    options->add("-snes_monitor");

    options->add("-ts_monitor");
    options->add("-ts_error_if_step_fails");

} // addMonitoring


// ------------------------------------------------------------------------------------------------
// Add default solver tolerances to options.
void
pylith::utils::_PetscOptions::addSolverTolerances(PetscOptions* options) {
    assert(options);

    options->add("-ksp_rtol", "1.0e-12");
    options->add("-ksp_atol", "1.0e-12");
    options->add("-ksp_error_if_not_converged");

    options->add("-snes_rtol", "1.0e-12");
    options->add("-snes_atol", "1.0e-9");
    options->add("-snes_error_if_not_converged");

} // addSolverTolerances


// End of file
