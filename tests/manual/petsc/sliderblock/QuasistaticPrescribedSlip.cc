/**
 * Quasistatic simulation with prescribed slip.
 *
 * Fully implicit solve: F(t,s) = 0.
 *
 * solution = [u0, u1, u2, u3, lambda]
 *
 * -2*ka*u0 + ka*u1 = 0
 * ka*u0 - ka*u1 + lambda = 0
 * -kb*u2 + kb*u3 - lambda = 0
 * kb*u2 - 2*kb*u3 + kb*u4 = 0
 * u1 - u2 + d = 0
 */

#include <portinfo>

#include "QuasistaticPrescribedSlip.hh"

#include "PrescribedSlip.hh" // USES PrescribedSlip
#include "DisplacementBC.hh" // USES DisplacementBC

#include "petsc.h"

#include <cassert> // USES assert()

// --------------------------------------------------------------------------------------------------
QuasistaticPrescribedSlip::QuasistaticPrescribedSlip(void) {
    _hasLHSResidual = true;
    _hasLHSJacobian = true;
}


// --------------------------------------------------------------------------------------------------
QuasistaticPrescribedSlip::~QuasistaticPrescribedSlip(void) {}


// --------------------------------------------------------------------------------------------------
void
QuasistaticPrescribedSlip::_computeLHSResidual(const PetscReal t,
                                               const PetscVec solution,
                                               const PetscVec solutionDot,
                                               PetscVec residual) {
    PetscErrorCode err = 0;
    const PetscScalar* solutionArray = NULL;
    PetscScalar* residualArray = NULL;

    err = VecGetArrayRead(solution, &solutionArray);CHECK_ERROR(err);
    err = VecGetArray(residual,&residualArray);CHECK_ERROR(err);

    const PetscScalar* u = &solutionArray[0];
    const PetscScalar lambda = solutionArray[_numDOFAll-1];

    const PetscScalar d = SlipFnRamp::slip(t);
    const PetscScalar u4 = DisplacementBC::displacement(t);

    residualArray[0] = -2*_ka*u[0] + _ka*u[1];
    residualArray[1] = +_ka*u[0] - _ka*u[1] + lambda;
    residualArray[2] = -_kb*u[2] + _kb*u[3] - lambda;
    residualArray[3] = +_kb*u[2] - 2*_kb*u[3] + _kb*u4;
    residualArray[4] = u[1] - u[2] + d;

    err = VecRestoreArrayRead(solution, &solutionArray);CHECK_ERROR(err);
    err = VecRestoreArray(residual, &residualArray);CHECK_ERROR(err);
}


// --------------------------------------------------------------------------------------------------
void
QuasistaticPrescribedSlip::_computeLHSJacobian(const PetscReal t,
                                               const PetscVec solution,
                                               const PetscVec solutionDot,
                                               const PetscReal stshift,
                                               PetscMat jacobian,
                                               PetscMat preconditioner) {
    assert(5 == _numDOFAll);
    PetscInt indices[_numDOFAll];
    for (size_t i = 0; i < _numDOFAll; ++i) {
        indices[i] = i;
    } // for

    PetscScalar jacobianArray[_numDOFAll][_numDOFAll];
    for (size_t i = 0; i < _numDOFAll; ++i) {
        for (size_t j = 0; j < _numDOFAll; ++j) {
            jacobianArray[i][j] = 0.0;
        } // for
    } // for

    jacobianArray[0][0] = -2.0*_ka;
    jacobianArray[0][1] = +_ka;
    jacobianArray[1][0] = +_ka;
    jacobianArray[1][1] = -_ka;
    jacobianArray[1][4] = +1.0;
    jacobianArray[2][2] = -_kb;
    jacobianArray[2][3] = +_kb;
    jacobianArray[2][4] = -1.0;
    jacobianArray[3][2] = +_kb;
    jacobianArray[3][3] = -2.0*_kb;
    jacobianArray[4][1] = +1.0;
    jacobianArray[4][2] = -1.0;
    PetscErrorCode err = 0;
    err = MatSetValues(jacobian, _numDOFAll, indices, _numDOFAll,indices, &jacobianArray[0][0],
                       INSERT_VALUES);CHECK_ERROR(err);
}


// End of file
