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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/materials/RheologyPoroelasticity.hh" // implementation of object methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_DEBUG

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::RheologyPoroelasticity::RheologyPoroelasticity(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::RheologyPoroelasticity::~RheologyPoroelasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::RheologyPoroelasticity::deallocate(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::RheologyPoroelasticity::updateKernelConstants(pylith::real_array* kernelConstants,
                                                                 const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateKernelConstants(kernelConstants"<<kernelConstants<<", dt="<<dt<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables, implicit.
void
pylith::materials::RheologyPoroelasticity::addKernelsUpdateStateVarsImplicit(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                                                             const spatialdata::geocoords::CoordSys* coordsys,
                                                                             const bool _useStateVars) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVarsImplicit(kernels="<<kernels<<", coordsys="<<coordsys<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVarsImplicit


// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables, explicit.
void
pylith::materials::RheologyPoroelasticity::addKernelsUpdateStateVarsExplicit(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                                                             const spatialdata::geocoords::CoordSys* coordsys,
                                                                             const bool _useStateVars) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVarsExplicit(kernels="<<kernels<<", coordsys="<<coordsys<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVarsExplicit


// End of file
