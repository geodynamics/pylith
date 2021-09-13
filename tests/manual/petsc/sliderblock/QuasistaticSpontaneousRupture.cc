/**
 * Quasistatic simulation with spontaneous rupture.
 *
 * Fully implicit solve: F(t,s) = 0.
 *
 * solution = [u0, u1, u2, u3, v0, v1, v2, v3, lambda]
 *
 * -2*ka*u0 + ka*u1 = 0
 * ka*u0 - ka*u1 + (friction-lambda) = 0
 * -kb*u2 + kb*u3 - (friction-lambda) = 0
 * kb*u2 - 2*kb*u3 + kb*u4 = 0
 * \dot{u0} - v0 = 0
 * \dot{u1} - v1 = 0
 * \dot{u2} - v2 = 0
 * \dot{u3} - v3 = 0
 * lambda * (u2 - u1 - d) = 0
 */

#include <portinfo>

#include "QuasistaticSpontaneousRupture.hh"

#include "Friction.hh" // USES Friction
#include "DisplacementBC.hh" // USES DisplacementBC

#include "petsc.h"

#include <cassert> // USES assert()
#include <iostream> // USES std::cerr

// --------------------------------------------------------------------------------------------------
QuasistaticSpontaneousRupture::QuasistaticSpontaneousRupture(const char* friction) :
    _friction(NULL) {
    _isSpontaneous = true;
    _hasLHSResidual = true;
    _hasLHSJacobian = true;

    if (0 == strcasecmp(friction, "static")) {
	_friction = new StaticFriction();
    } else if (0 == strcasecmp(friction, "slip_weakening")) {
      _friction = new SlipWeakeningFriction(true);
    } else if (0 == strcasecmp(friction, "viscous")) {
      _friction = new ViscousFriction();
    } else if (0 == strcasecmp(friction, "rate_state")) {
      //_friction = new RateStateFriction(true);
    } else {
      std::cerr << "Unknown friction model '" << friction << "'. Using default (slip_weakening)." << std::endl;
      _friction = new SlipWeakeningFriction(true);      
    }
}


// --------------------------------------------------------------------------------------------------
QuasistaticSpontaneousRupture::~QuasistaticSpontaneousRupture(void) {}


// --------------------------------------------------------------------------------------------------
void
QuasistaticSpontaneousRupture::_setSolutionBounds(PetscTS ts) {
    PetscErrorCode err = 0;

    PetscVec lowerBound = NULL;
    err = VecCreate(PETSC_COMM_WORLD, &lowerBound);CHECK_ERROR(err);
    err = VecSetSizes(lowerBound, PETSC_DECIDE, _numDOFAll);CHECK_ERROR(err);
    err = VecSetFromOptions(lowerBound);CHECK_ERROR(err);
    err = VecSet(lowerBound, PETSC_NINFINITY);
    err = VecSetValue(lowerBound, _numDOFAll-1, 0.0, INSERT_VALUES);

    PetscVec upperBound = NULL;
    err = VecCreate(PETSC_COMM_WORLD, &upperBound);CHECK_ERROR(err);
    err = VecSetSizes(upperBound, PETSC_DECIDE, _numDOFAll);CHECK_ERROR(err);
    err = VecSetFromOptions(upperBound);CHECK_ERROR(err);
    err = VecSet(upperBound, PETSC_INFINITY);

    PetscSNES snes = NULL;
    err = TSGetSNES(ts, &snes);CHECK_ERROR(err);
    err = SNESSetType(snes, SNESVINEWTONRSLS);CHECK_ERROR(err);
    SNESLineSearch snesLineSearch = NULL;
    err = SNESGetLineSearch(snes, &snesLineSearch);CHECK_ERROR(err);
    err = SNESLineSearchSetType(snesLineSearch, SNESLINESEARCHBASIC);CHECK_ERROR(err);
    err = SNESSetFromOptions(snes);CHECK_ERROR(err);
    err = SNESVISetVariableBounds(snes, lowerBound, upperBound);CHECK_ERROR(err);

    err = VecSetValue(_solution, _numDOFAll-1, 1.0, INSERT_VALUES);
}


#include <iostream>
// --------------------------------------------------------------------------------------------------
void
QuasistaticSpontaneousRupture::_computeLHSResidual(const PetscReal t,
                                                   const PetscVec solution,
                                                   const PetscVec solutionDot,
                                                   PetscVec residual) {
  const PetscScalar zeroTolerance = 1.0e-12;
  
    PetscErrorCode err = 0;
    const PetscScalar* solutionArray = NULL;
    const PetscScalar* solutionDotArray = NULL;
    PetscScalar* residualArray = NULL;

    err = VecGetArrayRead(solution, &solutionArray);CHECK_ERROR(err);
    err = VecGetArrayRead(solutionDot, &solutionDotArray);CHECK_ERROR(err);
    err = VecGetArray(residual,&residualArray);CHECK_ERROR(err);

    const PetscScalar* u = &solutionArray[0];
    const PetscScalar* uDot = &solutionDotArray[0];
    const PetscScalar* v = &solutionArray[_numDOFDisp];
    const PetscScalar lambda = solutionArray[_numDOFAll-1];

    const PetscScalar u4 = DisplacementBC::displacement(t);

    const PetscScalar slip = u[2] - u[1];
    const PetscScalar slipRate = v[2] - v[1];
    const PetscScalar friction = _friction->traction(fabs(slip), fabs(slipRate));
    const PetscScalar d = _friction->lockedSlip();

    PetscScalar fc = friction - lambda;
    if (slipRate < -zeroTolerance) { fc *= -1; }

#if 1
    std::cout << "t:" << t
              << ", u1:" << u[1]
              << ", u2:" << u[2]
              << ", f:" << friction
              << ", l:" << lambda
              << ", slip:" << slip
              << ", slipRate:" << slipRate
              << ", fc:" << fc
              << std::endl;
#endif

    residualArray[0] = -2*_ka*u[0] + _ka*u[1];
    residualArray[1] = +_ka*u[0] - _ka*u[1] + fc;
    residualArray[2] = -_kb*u[2] + _kb*u[3] - fc;
    residualArray[3] = +_kb*u[2] - 2*_kb*u[3] + _kb*u4;

    residualArray[4] = uDot[0] - v[0];
    residualArray[5] = uDot[1] - v[1];
    residualArray[6] = uDot[2] - v[2];
    residualArray[7] = uDot[3] - v[3];
    
    residualArray[8] = lambda * (u[2] - u[1] - d);

    err = VecRestoreArrayRead(solution, &solutionArray);CHECK_ERROR(err);
    err = VecRestoreArrayRead(solutionDot, &solutionDotArray);CHECK_ERROR(err);
    err = VecRestoreArray(residual, &residualArray);CHECK_ERROR(err);
}


// --------------------------------------------------------------------------------------------------
void
QuasistaticSpontaneousRupture::_computeLHSJacobian(const PetscReal t,
                                                   const PetscVec solution,
                                                   const PetscVec solutionDot,
                                                   const PetscReal stshift,
                                                   PetscMat jacobian,
                                                   PetscMat preconditioner) {
    assert(9 == _numDOFAll);
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

    PetscErrorCode err = 0;
    const PetscScalar* solutionArray = NULL;
    const PetscScalar* solutionDotArray = NULL;

    err = VecGetArrayRead(solution, &solutionArray);CHECK_ERROR(err);
    err = VecGetArrayRead(solutionDot, &solutionDotArray);CHECK_ERROR(err);

    const PetscScalar* u = &solutionArray[0];
    const PetscScalar* v = &solutionArray[_numDOFDisp];
    const PetscScalar slip = u[2] - u[1];
    const PetscScalar slipRate = v[2] - v[1];
    const PetscScalar df = _friction->jacobianSlip(fabs(slip));
    const PetscScalar ddf = _friction->jacobianSlipRate(fabs(slipRate));
    const PetscScalar d = _friction->lockedSlip();

    const PetscScalar lambda = solutionArray[_numDOFAll-1];

    jacobianArray[0][0] = -2.0*_ka;
    jacobianArray[0][1] = +_ka;
    jacobianArray[1][0] = +_ka;
    jacobianArray[1][1] = -_ka + (lambda > 0.0 ? 0.0 : -df);
    jacobianArray[1][2] = (lambda > 0.0 ? 0.0 : +df);
    jacobianArray[1][5] = (lambda > 0.0 ? 0.0 : -ddf);
    jacobianArray[1][6] = (lambda > 0.0 ? 0.0 : +ddf);
    jacobianArray[1][8] = -1.0;
    jacobianArray[2][1] = (lambda > 0.0 ? 0.0 : +df);
    jacobianArray[2][2] = -_kb + (lambda > 0.0 ? 0.0 : -df);
    jacobianArray[2][3] = +_kb;
    jacobianArray[2][5] = (lambda > 0.0 ? 0.0 : +ddf);
    jacobianArray[2][6] = (lambda > 0.0 ? 0.0 : -ddf);
    jacobianArray[2][8] = +1.0;
    jacobianArray[3][2] = +_kb;
    jacobianArray[3][3] = -2.0*_kb;

    jacobianArray[4][0] = stshift;
    jacobianArray[4][4] = -1.0;
    jacobianArray[5][1] = stshift;
    jacobianArray[5][5] = -1.0;
    jacobianArray[6][2] = stshift;
    jacobianArray[6][6] = -1.0;
    jacobianArray[7][3] = stshift;
    jacobianArray[7][7] = -1.0;
    
    jacobianArray[8][1] = -lambda;
    jacobianArray[8][2] = +lambda;
    jacobianArray[8][8] = lambda > 0 ? u[2] - u[1] - d : 1.0;
    
    err = VecRestoreArrayRead(solution, &solutionArray);CHECK_ERROR(err);
    err = VecRestoreArrayRead(solutionDot, &solutionDotArray);CHECK_ERROR(err);

    err = MatSetValues(jacobian, _numDOFAll, indices, _numDOFAll,indices, &jacobianArray[0][0],
                       INSERT_VALUES);CHECK_ERROR(err);
}


// --------------------------------------------------------------------------------------------------
void
QuasistaticSpontaneousRupture::_updateState(void) {
    PetscErrorCode err = 0;
    const PetscScalar* solutionArray = NULL;

    err = VecGetArrayRead(_solution, &solutionArray);CHECK_ERROR(err);

    const PetscScalar* u = &solutionArray[0];
    const PetscScalar* v = &solutionArray[_numDOFDisp];

    const PetscScalar slip = u[2] - u[1];
    const PetscScalar slipRate = v[2] - v[1];

    assert(_friction);
    _friction->updateState(fabs(slip), fabs(slipRate));

    // Set lambda to 1.0 as initial guess for next time step.
    err = VecSetValue(_solution, _numDOFAll-1, 1.0, INSERT_VALUES);
    
    err = VecRestoreArrayRead(_solution, &solutionArray);CHECK_ERROR(err);
}

// End of file
