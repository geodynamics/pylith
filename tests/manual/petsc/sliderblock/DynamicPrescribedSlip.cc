/**
 * Dynamic simulation with prescribed slip.
 *
 * Index 2 DAE solved using implicit-explicit (IMEX) method.
 *
 * Displacement and velocity equations are explicit.
 * Fault slip constraint is implicit.
 *
 * solution = [u0, u1, u2, u3, v0, v1, v2, v3, lambda]
 *
 * \dot{u0} = v0
 * \dot{u1} = v1
 * \dot{u2} = v2
 * \dot{u3} = v3
 * \dot{v0} = 1/ma * (-2*ka*u0 + ka*u1)
 * \dot{v1} = 1/ma * (ka*u0 - ka*u1 + lambda)
 * \dot{v2} = 1/mb * (-kb*u2 + kb*u3 - lambda)
 * \dot{v3} = 1/mb * (kb*u2 - 2*kb*u3 + kb*u4)
 * 1/ma * (ka*u0 - ka*u1) + 1/mb * (kb*u2 - kb*u3) + (1/ma + 1/mb) * lambda + \ddot{d} = 0
 */

#include <portinfo>

#include "DynamicPrescribedSlip.hh"

#include "PrescribedSlip.hh" // USES PrescribedSlip
#include "DisplacementBC.hh" // USES DisplacementBC

#include "petsc.h"

#include <cassert> // USES assert()

// --------------------------------------------------------------------------------------------------
DynamicPrescribedSlip::DynamicPrescribedSlip(void) {
    _hasRHSResidual = true;
    _hasLHSResidual = true;
    _hasLHSJacobian = true;
    _isDynamic = true;
    _tsAlgorithm = TSARKIMEX;
}


// --------------------------------------------------------------------------------------------------
DynamicPrescribedSlip::~DynamicPrescribedSlip(void) {}


// --------------------------------------------------------------------------------------------------
void
DynamicPrescribedSlip::_computeLHSResidual(const PetscReal t,
                                           const PetscVec solution,
                                           const PetscVec solutionDot,
                                           PetscVec residual) {
    PetscErrorCode err = 0;
    const PetscScalar* solutionArray = NULL;
    const PetscScalar* solutionDotArray = NULL;
    PetscScalar* residualArray = NULL;

    err = VecGetArrayRead(solution, &solutionArray);CHECK_ERROR(err);
    err = VecGetArrayRead(solutionDot, &solutionDotArray);CHECK_ERROR(err);
    err = VecGetArray(residual,&residualArray);CHECK_ERROR(err);

    const PetscScalar* u = &solutionArray[0];
    const PetscScalar lambda = solutionArray[_numDOFAll-1];
    const PetscScalar* udot = &solutionDotArray[0];
    const PetscScalar* vdot = &solutionDotArray[_numDOFDisp];

    const PetscScalar d2 = PrescribedSlip::slipAcc(t);

    for (size_t i = 0; i < _numDOFDisp; ++i) {
        residualArray[i] = udot[i];
    }
    for (size_t i = 0; i < _numDOFDisp; ++i) {
        residualArray[_numDOFDisp+i] = vdot[i];
    }

    residualArray[_numDOFAll-1] = (1.0/_ma + 1.0/_mb) * lambda
                                  + 1.0/_ma * (_ka*u[0] - _ka*u[1])
                                  + 1.0/_mb * (_kb*u[2] - _kb*u[3])
                                  + d2;

    const PetscScalar d = PrescribedSlip::slip(t);
    residualArray[_numDOFAll-1] += u[1] - u[2] + d;

    err = VecRestoreArrayRead(solution, &solutionArray);CHECK_ERROR(err);
    err = VecRestoreArrayRead(solutionDot, &solutionDotArray);CHECK_ERROR(err);
    err = VecRestoreArray(residual, &residualArray);CHECK_ERROR(err);
}


// --------------------------------------------------------------------------------------------------
void
DynamicPrescribedSlip::_computeRHSResidual(const PetscReal t,
                                           const PetscVec solution,
                                           PetscVec residual) {
    PetscErrorCode err = 0;
    const PetscScalar* solutionArray = NULL;
    PetscScalar* residualArray = NULL;

    err = VecGetArrayRead(solution, &solutionArray);CHECK_ERROR(err);
    err = VecGetArray(residual, &residualArray);CHECK_ERROR(err);

    const PetscScalar* u = &solutionArray[0];
    const PetscScalar* v = &solutionArray[_numDOFDisp];
    const PetscScalar lambda = solutionArray[_numDOFAll-1];

    const PetscScalar u4 = DisplacementBC::displacement(t);

    for (size_t i = 0; i < _numDOFDisp; ++i) {
        residualArray[i] = v[i];
    }

    residualArray[4] = 1.0/_ma*(-2*_ka*u[0] + _ka*u[1]);
    residualArray[5] = 1.0/_ma*(+_ka*u[0] - _ka*u[1] + lambda);
    residualArray[6] = 1.0/_mb*(-_kb*u[2] + _kb*u[3] - lambda);
    residualArray[7] = 1.0/_mb*(+_kb*u[2] - 2*_kb*u[3] + _kb*u4);

    residualArray[_numDOFAll-1] = 0.0;

    err = VecRestoreArrayRead(solution, &solutionArray);CHECK_ERROR(err);
    err = VecRestoreArray(residual, &residualArray);CHECK_ERROR(err);
}


// --------------------------------------------------------------------------------------------------
void
DynamicPrescribedSlip::_computeLHSJacobian(const PetscReal t,
                                           const PetscVec solution,
                                           const PetscVec solutionDot,
                                           const PetscReal shift,
                                           PetscMat jacobian,
                                           PetscMat preconditioner) {
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

    for (size_t i = 0; i < 2*_numDOFDisp; ++i) {
        jacobianArray[i][i] = shift;
    } // for

    jacobianArray[8][0] = +_ka/_ma;
    jacobianArray[8][1] = -_ka/_ma;
    jacobianArray[8][2] = +_kb/_mb;
    jacobianArray[8][3] = -_kb/_mb;
    jacobianArray[8][8] = 1.0/_ma + 1.0/_mb;

    jacobianArray[8][1] += +1.0;
    jacobianArray[8][2] += -1.0;

    PetscErrorCode err = 0;
    err = MatSetValues(jacobian, _numDOFAll, indices, _numDOFAll, indices, &jacobianArray[0][0],
                       INSERT_VALUES);CHECK_ERROR(err);
}


// End of file
