#include <portinfo>

#include "Formulation.hh"

#include "petscts.h"
#include "petscviewerhdf5.h"

#include <cassert> // USES assert()

//const double Formulation::_duration = 50.0;
//const double Formulation::_dtInitial = 0.01;
const double Formulation::_duration = 20.0;
const double Formulation::_dtInitial = 0.1;

const double Formulation::_ka = 1.0;
const double Formulation::_ma = 1.0;

const double Formulation::_kb = 1.0;
const double Formulation::_mb = 1.5;

const size_t Formulation::_numDOFDisp = 4;

// --------------------------------------------------------------------------------------------------
Formulation::Formulation(void) :
    _hasLHSResidual(false),
    _hasRHSResidual(false),
    _hasLHSJacobian(false),
    _hasRHSJacobian(false),
    _isDynamic(false),
    _isSpontaneous(false),
    _tsAlgorithm(TSBEULER),
    _ts(NULL),
    _tstamp(NULL),
    _viewer(NULL),
    _solution(NULL),
    _jacobianLHS(NULL),
    _jacobianRHS(NULL)
{}


// --------------------------------------------------------------------------------------------------
Formulation::~Formulation(void) {
    TSDestroy(&_ts);
    MatDestroy(&_jacobianLHS);
    MatDestroy(&_jacobianRHS);
    VecDestroy(&_solution);
    PetscViewerDestroy(&_viewer);
    VecDestroy(&_tstamp);
}


// --------------------------------------------------------------------------------------------------
void
Formulation::_createMatsAndVecs(void) {
    PetscErrorCode err = 0;
    err = MatDestroy(&_jacobianLHS);CHECK_ERROR(err);
    err = MatDestroy(&_jacobianRHS);CHECK_ERROR(err);
    err = VecDestroy(&_solution);CHECK_ERROR(err);

    _numDOFAll = (_isDynamic || _isSpontaneous) ? 2*_numDOFDisp + 1 : _numDOFDisp + 1;

    if (_hasLHSJacobian) {
        err = MatCreate(PETSC_COMM_WORLD, &_jacobianLHS);CHECK_ERROR(err);
        err = MatSetSizes(_jacobianLHS, PETSC_DECIDE, PETSC_DECIDE, _numDOFAll, _numDOFAll);CHECK_ERROR(err);
        err = MatSetFromOptions(_jacobianLHS);CHECK_ERROR(err);
        err = MatSetUp(_jacobianLHS);CHECK_ERROR(err);
    } // if

    if (_hasRHSJacobian) {
        err = MatCreate(PETSC_COMM_WORLD, &_jacobianRHS);CHECK_ERROR(err);
        err = MatSetSizes(_jacobianRHS, PETSC_DECIDE, PETSC_DECIDE, _numDOFAll, _numDOFAll);CHECK_ERROR(err);
        err = MatSetFromOptions(_jacobianRHS);CHECK_ERROR(err);
        err = MatSetUp(_jacobianRHS);CHECK_ERROR(err);
    } // if

    err = VecCreate(PETSC_COMM_WORLD, &_solution);CHECK_ERROR(err);
    err = VecSetSizes(_solution, PETSC_DECIDE, _numDOFAll);CHECK_ERROR(err);
    err = VecSetFromOptions(_solution);CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)_solution, "solution");CHECK_ERROR(err);
}


// --------------------------------------------------------------------------------------------------
void
Formulation::initialize(const char* const outputFilename) {
    _createMatsAndVecs();

    PetscErrorCode err = 0;
    err = TSDestroy(&_ts);CHECK_ERROR(err);
    err = TSCreate(PETSC_COMM_WORLD, &_ts);CHECK_ERROR(err);assert(_ts);
    err = TSSetType(_ts, _tsAlgorithm);CHECK_ERROR(err);
    err = TSSetExactFinalTime(_ts, TS_EXACTFINALTIME_STEPOVER);CHECK_ERROR(err); // Ok to step over final time.
    err = TSSetFromOptions(_ts);CHECK_ERROR(err);
    err = TSSetApplicationContext(_ts, (void*)this);CHECK_ERROR(err);
    err = TSSetProblemType(_ts, TS_NONLINEAR);CHECK_ERROR(err);
    err = TSSetTime(_ts, 0.0);CHECK_ERROR(err);
    err = TSSetTimeStep(_ts, _dtInitial);CHECK_ERROR(err);
    err = TSSetMaxTime(_ts, _duration);CHECK_ERROR(err);

    err = TSSetSolution(_ts, _solution);CHECK_ERROR(err);
    err = TSSetPostStep(_ts, poststep);CHECK_ERROR(err);

    if (_hasLHSResidual) { err = TSSetIFunction(_ts, NULL, computeLHSResidual, (void*)this);CHECK_ERROR(err); }
    if (_hasRHSResidual) { err = TSSetRHSFunction(_ts, NULL, computeRHSResidual, (void*)this);CHECK_ERROR(err); }
    if (_hasLHSJacobian) { err = TSSetIJacobian(_ts, _jacobianLHS, _jacobianLHS, computeLHSJacobian, (void*)this);CHECK_ERROR(err);}
    if (_hasRHSJacobian) { err = TSSetRHSJacobian(_ts, _jacobianRHS, _jacobianRHS, computeRHSJacobian, (void*)this);CHECK_ERROR(err);}

    _setSolutionBounds(_ts);

    const int tsize = 1;
    err = VecCreateMPI(PETSC_COMM_WORLD, tsize, 1, &_tstamp);CHECK_ERROR(err);assert(_tstamp);
    err = VecSetBlockSize(_tstamp, 1);CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) _tstamp, "time");CHECK_ERROR(err);
    err = PetscViewerHDF5Open(PETSC_COMM_WORLD, outputFilename, FILE_MODE_WRITE, &_viewer);CHECK_ERROR(err);
    err = PetscViewerHDF5SetBaseDimension2(_viewer, PETSC_TRUE);CHECK_ERROR(err);
}


// --------------------------------------------------------------------------------------------------
void
Formulation::solve(void) {
    PetscErrorCode err = TSSolve(_ts, NULL);CHECK_ERROR(err);
}


// --------------------------------------------------------------------------------------------------
void
Formulation::_poststep(void) {
    PetscFunctionBeginUser;

    this->_updateState();
    
    PetscErrorCode err = 0;
    PetscReal t = 0.0;
    PetscInt tindex = 0;
    err = TSGetTime(_ts, &t);CHECK_ERROR(err);
    err = TSGetStepNumber(_ts, &tindex);CHECK_ERROR(err);
 
    // Set time stamp
    err = VecSetValue(_tstamp, 0, t, INSERT_VALUES);CHECK_ERROR(err);
    err = VecAssemblyBegin(_tstamp);CHECK_ERROR(err);
    err = VecAssemblyEnd(_tstamp);CHECK_ERROR(err);

    // Write time stamp and solution
    err = PetscViewerHDF5PushGroup(_viewer, "/");CHECK_ERROR(err);
    err = PetscViewerHDF5PushTimestepping(_viewer);CHECK_ERROR(err);
    err = PetscViewerHDF5SetTimestep(_viewer, tindex);CHECK_ERROR(err);
    err = VecView(_tstamp, _viewer);CHECK_ERROR(err);
    err = VecView(_solution, _viewer);CHECK_ERROR(err);
    err = PetscViewerHDF5PopTimestepping(_viewer);CHECK_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer);CHECK_ERROR(err);

    PetscFunctionReturnVoid();
}


// --------------------------------------------------------------------------------------------------
PetscErrorCode
Formulation::computeLHSResidual(PetscTS ts,
                                PetscReal t,
                                PetscVec solution,
                                PetscVec solutionDot,
                                PetscVec residual,
                                void *context) {
    PetscFunctionBeginUser;

    Formulation* formulation = (Formulation*)context;assert(formulation);
    assert(formulation->_hasLHSResidual);
    formulation->_computeLHSResidual(t, solution, solutionDot, residual);

    PetscFunctionReturn(0);
}


// --------------------------------------------------------------------------------------------------
PetscErrorCode
Formulation::computeRHSResidual(PetscTS ts,
                                PetscReal t,
                                PetscVec solution,
                                PetscVec residual,
                                void *context) {
    PetscFunctionBeginUser;

    Formulation* formulation = (Formulation*)context;assert(formulation);
    assert(formulation->_hasRHSResidual);
    formulation->_computeRHSResidual(t, solution, residual);

    PetscFunctionReturn(0);
}


// --------------------------------------------------------------------------------------------------
PetscErrorCode
Formulation::computeLHSJacobian(PetscTS ts,
                                PetscReal t,
                                PetscVec solution,
                                PetscVec solutionDot,
                                PetscReal shift,
                                PetscMat jacobian,
                                PetscMat preconditioner,
                                void *context) {
    PetscFunctionBeginUser;

    Formulation* formulation = (Formulation*)context;assert(formulation);
    assert(formulation->_hasLHSJacobian);
    formulation->_computeLHSJacobian(t, solution, solutionDot, shift, jacobian, preconditioner);

    PetscErrorCode err = 0;
    err = MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY);CHECK_ERROR(err);
    err = MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY);CHECK_ERROR(err);
    if (jacobian != preconditioner) {
        err = MatAssemblyBegin(preconditioner, MAT_FINAL_ASSEMBLY);CHECK_ERROR(err);
        err = MatAssemblyEnd(preconditioner, MAT_FINAL_ASSEMBLY);CHECK_ERROR(err);
    } // if

    PetscFunctionReturn(0);
}


// --------------------------------------------------------------------------------------------------
PetscErrorCode
Formulation::computeRHSJacobian(PetscTS ts,
                                PetscReal t,
                                PetscVec solution,
                                PetscMat jacobian,
                                PetscMat preconditioner,
                                void *context) {
    PetscFunctionBeginUser;

    Formulation* formulation = (Formulation*)context;assert(formulation);
    assert(formulation->_hasRHSJacobian);
    formulation->_computeRHSJacobian(t, solution, jacobian, preconditioner);

    PetscErrorCode err = 0;
    err = MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY);CHECK_ERROR(err);
    err = MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY);CHECK_ERROR(err);
    if (jacobian != preconditioner) {
        err = MatAssemblyBegin(preconditioner, MAT_FINAL_ASSEMBLY);CHECK_ERROR(err);
        err = MatAssemblyEnd(preconditioner, MAT_FINAL_ASSEMBLY);CHECK_ERROR(err);
    } // if

    PetscFunctionReturn(0);
}


// --------------------------------------------------------------------------------------------------
PetscErrorCode
Formulation::poststep(PetscTS ts) {
    PetscFunctionBeginUser;

    Formulation* formulation = NULL;
    PetscErrorCode err = TSGetApplicationContext(ts, (void*)&formulation);CHECK_ERROR(err);assert(formulation);
    formulation->_poststep();

    PetscFunctionReturn(0);
}


// --------------------------------------------------------------------------------------------------
void
Formulation::_computeLHSResidual(const PetscReal,
                                 const PetscVec,
                                 const PetscVec,
                                 PetscVec) {
    throw std::logic_error("_computeLHSResidual() not implemented.");
}


// --------------------------------------------------------------------------------------------------
void
Formulation::_computeLHSJacobian(const PetscReal,
                                 const PetscVec,
                                 const PetscVec,
                                 const PetscReal,
                                 PetscMat,
                                 PetscMat) {
    throw std::logic_error("_computeLHSJacobian() not implemented.");
}


// --------------------------------------------------------------------------------------------------
void
Formulation::_computeRHSResidual(const PetscReal,
                                 const PetscVec,
                                 PetscVec) {
    throw std::logic_error("_computeRHSResidual() not implemented.");
}


// --------------------------------------------------------------------------------------------------
void
Formulation::_computeRHSJacobian(const PetscReal,
                                 const PetscVec,
                                 PetscMat,
                                 PetscMat) {
    throw std::logic_error("_computeRHSJacobian() not implemented.");
}


// --------------------------------------------------------------------------------------------------
void
Formulation::_updateState(void) {}

// --------------------------------------------------------------------------------------------------
void
Formulation::_setSolutionBounds(PetscTS ts) {}


// End of file
