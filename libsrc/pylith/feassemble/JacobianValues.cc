// -*- C++ -*-
//
// ------------------------------------------------------------------------------------------------
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
// ------------------------------------------------------------------------------------------------
//

#include <portinfo>

#include "JacobianValues.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/DSLabelAccess.hh" // USES DSLabelAccess

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

namespace pylith {
    namespace feassemble {
        class _JacobianValues {
public:

            static
            void evaluateKernel(scalar_array* cellMat,
                                const JacobianValues::JacobianKernel& kernel,
                                const PylithReal t,
                                const PylithReal dt,
                                const PylithReal s_tshift,
                                const pylith::topology::Field& solution,
                                const pylith::feassemble::DSLabelAccess& dsLabel);

        };
    }
}

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::JacobianValues::JacobianValues(void) {
    GenericComponent::setName("jacobianvalues");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::JacobianValues::~JacobianValues(void) {}


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::JacobianValues::setKernels(const std::vector<JacobianKernel>& kernelsJacobian,
                                               const std::vector<JacobianKernel>& kernelsPrecond) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernels(# Jacobian kernels="<<kernelsJacobian.size()<<", # preconditioner kernels="<<kernelsPrecond.size()<<")");

    _kernelsJacobian = kernelsJacobian;
    _kernelsPrecond = kernelsPrecond;

    PYLITH_METHOD_END;
} // setKernels


// ------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::JacobianValues::computeLHSJacobian(PetscMat jacobianMat,
                                                       PetscMat precondMat,
                                                       const PylithReal t,
                                                       const PylithReal dt,
                                                       const PylithReal s_tshift,
                                                       const pylith::topology::Field& solution,
                                                       const pylith::feassemble::DSLabelAccess& dsLabel) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSJacobian (jacobianMat = "<<jacobianMat<<", precondMat = "<<precondMat<<", t = "<<t<<", dt = "<<dt<<", solution = "<<solution.getLabel()<<") ");

    PetscDM dm = dsLabel.dm();assert(dm);

    PetscErrorCode err;
    const PetscInt numCells = dsLabel.numCells();
    const PetscInt* cellIndices = NULL;
    err = ISGetIndices(dsLabel.cellsIS(), &cellIndices);PYLITH_CHECK_ERROR(err);

    PetscInt totalDof = 0;
    err = PetscDSGetTotalDimension(dsLabel.ds(), &totalDof);PYLITH_CHECK_ERROR(err);
    scalar_array cellMat(totalDof*totalDof);

    if (jacobianMat) {
        for (PetscInt iCell = 0; iCell < numCells; ++iCell) {
            const PetscInt cell = cellIndices[iCell];
            cellMat = 0.0;

            for (size_t i = 0; i < _kernelsJacobian.size(); ++i) {
                _JacobianValues::evaluateKernel(&cellMat, _kernelsJacobian[i], t, dt, s_tshift, solution, dsLabel);
            } // for
            err = DMPlexMatSetClosure(dsLabel.dm(), NULL, NULL, jacobianMat, cell, &cellMat[0],
                                      INSERT_VALUES);PYLITH_CHECK_ERROR(err);
        } // for
    } // if

    if (precondMat && (jacobianMat != precondMat)) {
        for (PetscInt iCell = 0; iCell < numCells; ++iCell) {
            const PetscInt cell = cellIndices[iCell];
            cellMat = 0.0;

            for (size_t i = 0; i < _kernelsPrecond.size(); ++i) {
                _JacobianValues::evaluateKernel(&cellMat, _kernelsPrecond[i], t, dt, s_tshift, solution, dsLabel);
            } // for
            err = DMPlexMatSetClosure(dsLabel.dm(), NULL, NULL, precondMat, cell, &cellMat[0],
                                      INSERT_VALUES);PYLITH_CHECK_ERROR(err);
        } // for
    } // if
    err = ISRestoreIndices(dsLabel.cellsIS(), &cellIndices);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ------------------------------------------------------------------------------------------------
// Compute Jacobian values associated with s_tshift on diagonal.
void
pylith::feassemble::JacobianValues::blockDiag_tshift(pylith::scalar_array* cellMat,
                                                     const PylithReal t,
                                                     const PylithReal dt,
                                                     const PylithReal s_tshift,
                                                     const PylithInt trialDof,
                                                     const PylithInt trialOff,
                                                     const PylithInt basisDof,
                                                     const PylithInt basisOff,
                                                     const PylithInt totalDim) {
    assert(trialDof == basisDof);

    const PylithScalar value = s_tshift;
    for (PetscInt iDof = 0; iDof < trialDof; ++iDof) {
        const PylithInt index = (trialOff+iDof) * totalDim + basisOff + iDof;
        (*cellMat)[index] = value;
    } // for
} // blockDiag_tshift


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::_JacobianValues::evaluateKernel(scalar_array* cellMat,
                                                    const JacobianValues::JacobianKernel& kernel,
                                                    const PylithReal t,
                                                    const PylithReal dt,
                                                    const PylithReal s_tshift,
                                                    const pylith::topology::Field& solution,
                                                    const pylith::feassemble::DSLabelAccess& dsLabel) {
    const pylith::topology::Field::SubfieldInfo& infoTrial = solution.getSubfieldInfo(kernel.subfieldTrial.c_str());
    const pylith::topology::Field::SubfieldInfo& infoBasis = solution.getSubfieldInfo(kernel.subfieldBasis.c_str());
    const size_t i_trial = infoTrial.index;
    const size_t i_basis = infoBasis.index;

    PetscErrorCode err;
    PetscInt trialOff, trialDof, basisOff, basisDof;
    PetscFE fe = NULL;
    err = PetscDSGetFieldOffset(dsLabel.ds(), i_trial, &trialOff);PYLITH_CHECK_ERROR(err);
    err = PetscDSGetDiscretization(dsLabel.ds(), i_trial, (PetscObject*) &fe);PYLITH_CHECK_ERROR(err);
    err = PetscFEGetDimension(fe, &trialDof);PYLITH_CHECK_ERROR(err);

    err = PetscDSGetFieldOffset(dsLabel.ds(), i_basis, &basisOff);PYLITH_CHECK_ERROR(err);
    err = PetscDSGetDiscretization(dsLabel.ds(), i_basis, (PetscObject*) &fe);PYLITH_CHECK_ERROR(err);
    err = PetscFEGetDimension(fe, &basisDof);PYLITH_CHECK_ERROR(err);

    PetscInt totalDof = 0;
    err = PetscDSGetTotalDimension(dsLabel.ds(), &totalDof);PYLITH_CHECK_ERROR(err);

    kernel.function(cellMat, t, dt, s_tshift, trialDof, trialOff, basisDof, basisOff, totalDof);
} // evaluateKernel


// End of file
