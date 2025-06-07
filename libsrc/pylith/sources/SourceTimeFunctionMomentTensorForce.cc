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

#include "pylith/sources/SourceTimeFunctionMomentTensorForce.hh" // implementation of object methods

#include "pylith/feassemble/Integrator.hh" // USES NEW_JACOBIAN_NEVER

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_DEBUG

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::SourceTimeFunctionMomentTensorForce::SourceTimeFunctionMomentTensorForce(void) :
    _JacobianTriggers(pylith::feassemble::Integrator::NEW_JACOBIAN_NEVER) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::SourceTimeFunctionMomentTensorForce::~SourceTimeFunctionMomentTensorForce(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::sources::SourceTimeFunctionMomentTensorForce::deallocate(void) {
}


// ---------------------------------------------------------------------------------------------------------------------
// Get triggers for needing to compute the elastic constants for the RHS Jacobian.
int
pylith::sources::SourceTimeFunctionMomentTensorForce::getJacobianTriggers(void) const {
    return _JacobianTriggers;
} // getJacobianTriggers


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::sources::SourceTimeFunctionMomentTensorForce::updateKernelConstants(pylith::real_array* kernelConstants,
                                                                            const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateKernelConstants(kernelConstants"<<kernelConstants<<", dt="<<dt<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables.
void
pylith::sources::SourceTimeFunctionMomentTensorForce::addKernelsUpdateStateVars(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                                                                const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVars(kernels="<<kernels<<", coordsys="<<typeid(coordsys).name()<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVars


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary field values.
void
pylith::sources::SourceTimeFunctionMomentTensorForce::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                                           const PylithReal t,
                                                                           const PylithReal timeScale,
                                                                           spatialdata::units::Nondimensional* _normalizer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateAuxiliaryField(auxiliaryField"<<auxiliaryField<<", t="<<t<<", timeScale"<<timeScale<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // updateKernelConstants


// End of file
