// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, Rice University
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "IntegratorPointwise.hh" // implementation of class methods

#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> \
    // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorPointwise::IntegratorPointwise(void) :
    _normalizer(new spatialdata::units::Nondimensional),
    _gravityField(NULL),
    _auxField(NULL),
    _derivedField(NULL),
    _logger(NULL),
    _needNewRHSJacobian(true),
    _needNewLHSJacobian(true)
{}

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorPointwise::~IntegratorPointwise(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorPointwise::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    Observers::deallocate();

    delete _normalizer; _normalizer = NULL;
    delete _logger; _logger = NULL;
    delete _auxField; _auxField = NULL;
    delete _derivedField; _derivedField = NULL;

    _gravityField = NULL; // :KLUDGE: Memory managed by Python object. :TODO: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Get auxiliary field.
const pylith::topology::Field*
pylith::feassemble::IntegratorPointwise::auxField(void) const {
    return _auxField;
} // auxField

// ----------------------------------------------------------------------
// Get derived field.
const pylith::topology::Field*
pylith::feassemble::IntegratorPointwise::derivedField(void) const {
    return _derivedField;
} // derivedField

// ----------------------------------------------------------------------
// Set database for auxiliary fields.
void
pylith::feassemble::IntegratorPointwise::auxFieldDB(spatialdata::spatialdb::SpatialDB* value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("auxFieldDB(value="<<value<<")");

    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory(); assert(factory);
    factory->queryDB(value);

    PYLITH_METHOD_END;
} // auxFieldDB


// ----------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::feassemble::IntegratorPointwise::auxSubfieldDiscretization(const char* name,
                                                                   const int basisOrder,
                                                                   const int quadOrder,
                                                                   const bool isBasisContinuous,
                                                                   const pylith::topology::FieldBase::SpaceEnum feSpace) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("auxSubfieldDiscretization(name="<<name<<", basisOrder="<<basisOrder<<", quadOrder="<<quadOrder<<", isBasisContinuous="<<isBasisContinuous<<")");

    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory(); assert(factory);
    factory->subfieldDiscretization(name, basisOrder, quadOrder, isBasisContinuous, feSpace);

    PYLITH_METHOD_END;
} // auxSubfieldDiscretization


// ----------------------------------------------------------------------
// Check whether RHS Jacobian needs to be recomputed.
bool
pylith::feassemble::IntegratorPointwise::needNewRHSJacobian(void) const {
    return _needNewRHSJacobian;
} // needNewRHSJacobian

// ----------------------------------------------------------------------
// Check whether LHS Jacobian needs to be recomputed.
bool
pylith::feassemble::IntegratorPointwise::needNewLHSJacobian(void) const {
    return _needNewLHSJacobian;
} // needNewLHSJacobian

// ----------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::feassemble::IntegratorPointwise::normalizer(const spatialdata::units::Nondimensional& dim) {
    PYLITH_COMPONENT_DEBUG("normalizer(dim="<<typeid(dim).name()<<")");

    if (!_normalizer) {
        _normalizer = new spatialdata::units::Nondimensional(dim);
    } else {
        *_normalizer = dim;
    } // if/else
} // normalizer


// ----------------------------------------------------------------------
// Set gravity field.
void
pylith::feassemble::IntegratorPointwise::gravityField(spatialdata::spatialdb::GravityField* const g) {
    _gravityField = g;
} // gravityField


// ----------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::feassemble::IntegratorPointwise::prestep(const PylithReal t,
                                                 const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("prestep(t="<<t<<", dt="<<dt<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // prestep


// ----------------------------------------------------------------------
// Update auxiliary fields at end of time step.
void
pylith::feassemble::IntegratorPointwise::poststep(const PylithReal t,
                                                  const PylithInt tindex,
                                                  const PylithReal dt,
                                                  const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("poststep(t="<<t<<", dt="<<dt<<") empty method");

    _updateStateVars(t, dt, solution);

    const bool infoOnly = false;
    notifyObservers(t, tindex, solution, infoOnly);

    PYLITH_METHOD_END;
} // poststep


// ----------------------------------------------------------------------
// Set constants used in finite-element integrations.
void
pylith::feassemble::IntegratorPointwise::_setFEConstants(const pylith::topology::Field& solution,
                                                         const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEConstants(solution="<<solution.label()<<", dt="<<dt<<")");

    PetscDS prob = NULL;
    PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);

    // Pointwise functions have been set in DS
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err); assert(prob);
    err = PetscDSSetConstants(prob, 0, NULL); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEConstants


// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorPointwise::_updateStateVars(const PylithReal t,
                                                          const PylithReal dt,
                                                          const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateStateVars(t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    if (0 == _updateStateVarsKernels.size()) {
        PYLITH_METHOD_END;
    } // if

    assert(_auxField);

    PetscErrorCode err;
#if 0
    PetscDM  auxDM = _auxField->dmMesh(), stateVarDM, dms[2], superDM;
    PetscIS  stateVarIS, *superIS;
    PetscVec stateVarVec, A, locA;

    // HAPPEN ONCE
    // Create subDM holding only the state vars
    err = DMCreateSubDM(auxDM, numStateSubfields, stateSubfields, &stateVarIS, &stateVarDM);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(stateVarDM, &stateVarVec);PYLITH_CHECK_ERROR(err);
    // Create superDM of {state vars, solution}
    dms[0] = stateVarDM; dms[1] = solution->dmMesh();
    err = DMCreateSuperDM(dms, 2, &superIS, &superDM);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(superDM, &A);PYLITH_CHECK_ERROR(err);
    err = DMCreateLocalVector(superDM, &locA);PYLITH_CHECK_ERROR(err);
    
    // Copy state vars from auxiliary vars
    err = VecISCopy(stateVarVec, stateVarIS, SCATTER_FORWARD, _auxField->globalVector());PYLITH_CHECK_ERROR(err);
    // Copy current state vars and solution into superDM space
    err = VecISCopy(A, subis[0], SCATTER_FORWARD, stateVarVec);PYLITH_CHECK_ERROR(err);
    err = VecISCopy(A, subis[1], SCATTER_FORWARD, solution->globalVector());PYLITH_CHECK_ERROR(err);
    // Move superDM data to a local vector
    err = DMGlobalToLocalBegin(superDM, A, INSERT_VALUES, locA);PYLITH_CHECK_ERROR(err);
    err = DMGlobalToLocalEnd(superDM, A, INSERT_VALUES, locA);PYLITH_CHECK_ERROR(err);
    // Attach superDM data as auxiliary for updateStateVar kernel
    err = PetscObjectCompose((PetscObject) auxDM, "dmAux", (PetscObject) superDM);CHKERRQ(ierr);
    err = PetscObjectCompose((PetscObject) auxDM, "A",     (PetscObject) locA);CHKERRQ(ierr);
#else
    PetscDM dmState = _auxField->dmMesh();
    err = PetscObjectCompose((PetscObject) dmState, "dmAux", (PetscObject) solution.dmMesh());PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmState, "A", (PetscObject) solution.localVector());PYLITH_CHECK_ERROR(err);
#endif
    
    _setFEConstants(*_auxField, dt);

    // Set update kernel for each auxiliary subfield.
    const pylith::string_vector& subfieldNames = _auxField->subfieldNames();
    const size_t numSubfields = subfieldNames.size();
    PetscPointFunc* stateVarsKernels = (numSubfields > 0) ? new PetscPointFunc[numSubfields] : NULL;
    // By default, set all auxiliary subfield update kernels to NULL.
    for (size_t i = 0; i < numSubfields; ++i) {
        stateVarsKernels[i] = NULL;
    } // for
    for (UpdateStateVarsMap::const_iterator iter = _updateStateVarsKernels.begin(); iter != _updateStateVarsKernels.end(); ++iter) {
        const pylith::topology::Field::SubfieldInfo& sinfo = _auxField->subfieldInfo(iter->first.c_str());
        stateVarsKernels[sinfo.index] = iter->second;
    } // for

    err = DMProjectFieldLocal(dmState, t, _auxField->localVector(), stateVarsKernels, INSERT_VALUES, _auxField->localVector());PYLITH_CHECK_ERROR(err);
    delete[] stateVarsKernels; stateVarsKernels = NULL;

#if 0
    // HAPPEN ONCE
    // Destroy superDM stuff
    for (i = 0; i < 2; ++i) {err = ISDestroy(&superIS[i]);PYLITH_CHECK_ERROR(err);}
    err = DMDestroy(&superDM);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&stateVarIS);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&stateVarDM);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&locA);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&A);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&stateVarVec);PYLITH_CHECK_ERROR(err);
#endif
    
    PYLITH_METHOD_END;
} // _updateStateVars


// ----------------------------------------------------------------------
// Compute field derived from solution and auxiliary field.
void
pylith::feassemble::IntegratorPointwise::_computeDerivedFields(const PylithReal t,
                                                               const PylithReal dt,
                                                               const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_computeDerivedFields(t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");



    PYLITH_METHOD_END;
} // _computeDerivedFields

// End of file
