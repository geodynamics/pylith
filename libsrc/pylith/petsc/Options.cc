// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "Options.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/materials/Material.hh" // USES Material

#include "pylith/exceptions/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*
#include "pylith/utils/mpi.hh" // USES isRoot()

#include <cassert>

namespace pylith {
    namespace petsc {
        class _Options {
public:

            /// Write options to journal;
            static
            void write(pythia::journal::info_t& info,
                       const char* heading,
                       const Options& options);

            /** Check if simulation is running in parallel.
             *
             * @param[in] solution Solution field for problem.
             * @returns True if solving problem in parallel, False if in serial.
             */
            static
            bool isParallel(void);

            /** Add debugging options.
             *
             * @param[in] options PETSc options.
             */
            static
            void addTesting(Options* options);

            /** Add monitoring options.
             *
             * @param[in] options PETSc options.
             */
            static
            void addMonitoring(Options* options);

            /** Add collective I/O options.
             *
             * @param[in] options PETSc options.
             */
            static
            void addCollectiveIO(Options* options);

            /** Add default solver tolerances to options.
             *
             * @param[in] options PETSc options.
             */
            static
            void addSolverTolerances(Options* options);

            /** Add initial guess options.
             *
             * @param[in] options PETSc options.
             */
            static
            void addInitialGuess(Options* options);

            /** Add adaptive time stepping options.
             *
             * @param[in] options PETSc options.
             */
            static
            void addAdaptiveTimeStepping(Options* options);

        };
    }
}

// ------------------------------------------------------------------------------------------------
const int pylith::petsc::Defaults::NONE = 0x0;
const int pylith::petsc::Defaults::MONITORS = 0x1;
const int pylith::petsc::Defaults::SOLVER = 0x2;
const int pylith::petsc::Defaults::PARALLEL = 0x4;
const int pylith::petsc::Defaults::INITIAL_GUESS = 0x8;
const int pylith::petsc::Defaults::TESTING = 0x10;
const int pylith::petsc::Defaults::COLLECTIVE_IO = 0x20;
const int pylith::petsc::Defaults::TS_ADAPTIVE = 0x40;
const int pylith::petsc::Defaults::TS_IMPULSE = 0x80;

// ------------------------------------------------------------------------------------------------
// Set default PETSc solver options based on solution field and material.
void
pylith::petsc::Defaults::set(const pylith::materials::Material* material,
                             const bool hasFault,
                             const int flags) {
    PYLITH_METHOD_BEGIN;
    assert(material);

    if (!flags) {
        PYLITH_METHOD_END;
    } // if

    pylith::petsc::Options* options = NULL;
    if (flags & SOLVER) {
        const bool isParallel = flags & PARALLEL || _Options::isParallel();
        options = material->getSolverDefaults(isParallel, hasFault);
    } // if
    if (!options) {
        options = new pylith::petsc::Options();
    } // if
    assert(options);

    _Options::addSolverTolerances(options);
    if (flags & INITIAL_GUESS) {
        _Options::addInitialGuess(options);
    } // if
    if (flags & TESTING) {
        _Options::addTesting(options);
    } // if
    if (flags & MONITORS) {
        _Options::addMonitoring(options);
    } // if
    if (flags & COLLECTIVE_IO) {
        _Options::addCollectiveIO(options);
    } // if
    if (flags & TS_ADAPTIVE) {
        _Options::addAdaptiveTimeStepping(options);
    } // if

    options->set();
    delete options;options = NULL;

    PYLITH_METHOD_END;
} // set


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::petsc::Options::Options(void) {
    GenericComponent::setName("petscoptions");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::petsc::Options::~Options(void) {}


// ------------------------------------------------------------------------------------------------
// Add PETSc option.
void
pylith::petsc::Options::add(const char* name,
                            const char* value) {
    _options[std::string(name)] = value;
} // add


// ------------------------------------------------------------------------------------------------
// Remove PETSc option.
void
pylith::petsc::Options::remove(const char* name) {
    const options_t::iterator iter = _options.find(std::string(name));
    if (_options.end() != iter) {
        _options.erase(iter);
    } // if
} // remove


// ------------------------------------------------------------------------------------------------
// Clear PETSc options.
void
pylith::petsc::Options::clear(void) {
    _options.clear();
} // clear


// ------------------------------------------------------------------------------------------------
// Set PETSc options.
void
pylith::petsc::Options::set(void) {
    PYLITH_METHOD_BEGIN;

    pylith::petsc::Options optionsUsed;
    pylith::petsc::Options optionsIgnored;
    for (options_t::iterator iter = _options.begin(); iter != _options.end(); ++iter) {
        const char* name = iter->first.c_str();
        const char* value = iter->second.c_str();

        PetscBool exists = PETSC_FALSE;
        PylithCallPetsc(PetscOptionsHasName(NULL, NULL, name, &exists));
        if (!exists) {
            PylithCallPetsc(PetscOptionsSetValue(NULL, name, value));
            optionsUsed.add(name, value);
        } else {
            optionsIgnored.add(name, value);
        } // if/else
    } // for

    pythia::journal::info_t info(pylith::journal::application_flow);
    if (info.state() && pylith::utils::MPI::isRoot()) {
        _Options::write(info, "Setting PETSc options:", optionsUsed);
        if (optionsIgnored._options.size() > 0) {
            _Options::write(info, "Using user values rather then the following default PETSc options:", optionsIgnored);
        } // if
    } // if

    PYLITH_METHOD_END;
} // set


// ------------------------------------------------------------------------------------------------
// Set PETSc options, overriding any previously set options with the same name.
void
pylith::petsc::Options::override (void) {
    PYLITH_METHOD_BEGIN;

    for (options_t::iterator iter = _options.begin(); iter != _options.end(); ++iter) {
        const char* name = iter->first.c_str();
        const char* value = iter->second.c_str();

        PylithCallPetsc(PetscOptionsSetValue(NULL, name, value));
    } // for

    pythia::journal::info_t info(pylith::journal::application_flow);
    if (info.state()) {
        _Options::write(info, "Setting PETSc options:", (*this));
    } // if

    PYLITH_METHOD_END;
} // set

// ------------------------------------------------------------------------------------------------
// Write options to journal;
void
pylith::petsc::_Options::write(pythia::journal::info_t& info,
                               const char* heading,
                               const Options& options) {
    PYLITH_METHOD_BEGIN;

    info << pythia::journal::at(__HERE__)
         << heading << "\n";
    const Options::options_t::const_iterator begin = options._options.begin();
    const Options::options_t::const_iterator end = options._options.end();
    for (Options::options_t::const_iterator iter = begin; iter != end; ++iter) {
        const std::string name = iter->first.substr(1);
        const std::string value = iter->second;
        if (iter->second.empty()) {
            info << "    " << name << " = true\n";
        } else {
            info << "    " << name << " = " << value << "\n";
        } // if/else
    } // for
    info << pythia::journal::endl;

    PYLITH_METHOD_END;
} // write


// ------------------------------------------------------------------------------------------------
// Check if simulation is running in parallel.
bool
pylith::petsc::_Options::isParallel(void) {
    PYLITH_METHOD_BEGIN;

    MPI_Comm comm = PETSC_COMM_WORLD;
    int numProcs = 0;
    MPI_Comm_size(comm, &numProcs);

    PYLITH_METHOD_RETURN(numProcs > 1);
} // isParallel


// ------------------------------------------------------------------------------------------------
// Add debugging options.
void
pylith::petsc::_Options::addTesting(Options* options) {
    assert(options);

    // -checkstack only works with PetscInitialize()
    options->add("-malloc_dump");
} // setDebugging


// ------------------------------------------------------------------------------------------------
// Add monitoring options.
void
pylith::petsc::_Options::addMonitoring(Options* options) {
    assert(options);

    options->add("-ksp_error_if_not_converged");
    options->add("-snes_error_if_not_converged");
    options->add("-ts_error_if_step_fails");

    pythia::journal::info_t solver(pylith::journal::solver);
    solver.activate();
} // addMonitoring


// ------------------------------------------------------------------------------------------------
// Add collective I/O options.
void
pylith::petsc::_Options::addCollectiveIO(Options* options) {
    assert(options);

    options->add("-viewer_hdf5_collective");

} // addMonitoring


// ------------------------------------------------------------------------------------------------
// Add default solver tolerances to options.
void
pylith::petsc::_Options::addSolverTolerances(Options* options) {
    assert(options);

    options->add("-ksp_rtol", "1.0e-14");
    options->add("-ksp_atol", "1.0e-7");

    options->add("-snes_rtol", "1.0e-14");
    options->add("-snes_atol", "5.0e-7");
} // addSolverTolerances


// ------------------------------------------------------------------------------------------------
// Add initial guess defaults.
void
pylith::petsc::_Options::addInitialGuess(Options* options) {
    assert(options);

    options->add("-ksp_guess_type", "pod");
    options->add("-ksp_guess_pod_size", "8");

} // addInitialGuess


// ------------------------------------------------------------------------------------------------
// Add adaptive time stepping defaults.
void
pylith::petsc::_Options::addAdaptiveTimeStepping(Options* options) {
    assert(options);

    options->add("-ts_adapt_type", "basic");
    options->add("-ts_adapt_safety", "0.2");
    options->add("-ts_adapt_reject_safety", "0.1");
    options->add("-ts_atol", "0.05");
    options->add("-ts_rtol", "0.05");

    options->add("-ts_adapt_monitor", "true");


} // addAdaptiveTimeStepping


// End of file
