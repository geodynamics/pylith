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
#include "pylith/meshio/OutputManager.hh" // USES OutputManager

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> \
    // USES std::runtime_error

#include <iostream> // debugging

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorPointwise::IntegratorPointwise(void) :
    _normalizer(new spatialdata::units::Nondimensional),
    _gravityField(NULL),
    _auxField(NULL),
    _output(NULL),
    _logger(NULL),
    _needNewRHSJacobian(true),
    _needNewLHSJacobian(true)
{} // constructor

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

    delete _normalizer; _normalizer = NULL;
    delete _logger; _logger = NULL;
    delete _auxField; _auxField = NULL;

    _gravityField = NULL; // :KLUDGE: Memory managed by Python object. :TODO: Use shared pointer.
    _output = NULL; // :KLUDGE: Memory managed by Python object. :TODO: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Get auxiliary field.
const pylith::topology::Field&
pylith::feassemble::IntegratorPointwise::auxField(void) const
{
    PYLITH_METHOD_BEGIN;

    assert(_auxField);

    PYLITH_METHOD_RETURN(*_auxField);
} // auxField

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
// Set output manager.
void
pylith::feassemble::IntegratorPointwise::output(pylith::meshio::OutputManager* manager) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("output(manager="<<manager<<")");

    _output = manager;

    PYLITH_METHOD_END;
} // output


// ----------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::feassemble::IntegratorPointwise::prestep(const double t,
                                                 const double dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("prestep(t="<<t<<", dt="<<dt<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // prestep


// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorPointwise::updateStateVars(const PylithReal t,
                                                         const PylithReal dt,
                                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateStateVars(t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<") empty method");

    if (0 == _updateStateVarsKernels.size()) {
        PYLITH_METHOD_END;
    } // if

    assert(_auxField);

    PetscErrorCode err;
#if 1
    PetscDM  auxDM = _auxField->dmMesh(), solutionDM = solution.dmMesh(), stateVarDM, dms[2], superDM;
    PetscIS  stateVarIS, *superIS;
	//_auxField->view("AUXILIARY FIELD"); //:DEBUG:

    const pylith::string_vector& subfieldNames = _auxField->subfieldNames();
    const size_t numSubfields = subfieldNames.size();
	pylith::int_array stateSubfieldIndices(numSubfields);

	size_t numStateSubfields = 0;
	for (size_t iSubfield=0; iSubfield < numSubfields; ++iSubfield) {
		const pylith::topology::Field::SubfieldInfo& info = _auxField->subfieldInfo(subfieldNames[iSubfield].c_str());
		if (info.description.hasHistory) {
			stateSubfieldIndices[numStateSubfields++] = info.index;
		} // if
	} // for
	
    // HAPPEN ONCE
    // Create subDM holding only the state vars
    PetscVec stateVarVecGlobal = NULL;
    err = DMCreateSubDM(auxDM, numStateSubfields, &stateSubfieldIndices[0], &stateVarIS, &stateVarDM);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(stateVarDM, &stateVarVecGlobal);PYLITH_CHECK_ERROR(err);
	//err = DMView(stateVarDM, PETSC_VIEWER_STDOUT_WORLD); //:DEBUG:

    // Create superDM of {state vars, solution}
    PetscVec stateVarsSolutionGlobal = NULL, stateVarsSolutionLocal = NULL;
    dms[0] = stateVarDM; dms[1] = solutionDM;
    err = DMCreateSuperDM(dms, 2, &superIS, &superDM);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(superDM, &stateVarsSolutionGlobal);PYLITH_CHECK_ERROR(err);
    err = DMCreateLocalVector(superDM, &stateVarsSolutionLocal);PYLITH_CHECK_ERROR(err);
	//err = DMView(superDM, PETSC_VIEWER_STDOUT_WORLD); //:DEBUG:
    
    // Copy state vars from auxiliary vars
	PetscVec auxFieldVecGlobal = NULL;
	err = DMCreateGlobalVector(auxDM, &auxFieldVecGlobal);
    err = DMLocalToGlobalBegin(auxDM, _auxField->localVector(), INSERT_VALUES, auxFieldVecGlobal);PYLITH_CHECK_ERROR(err);
    err = DMLocalToGlobalEnd(auxDM, _auxField->localVector(), INSERT_VALUES, auxFieldVecGlobal);PYLITH_CHECK_ERROR(err);
    err = VecISCopy(auxFieldVecGlobal, stateVarIS, SCATTER_REVERSE, stateVarVecGlobal);PYLITH_CHECK_ERROR(err);
#if 1
	std::cout << std::endl << "stateVarVecGlobal" << std::endl; //:DEBUG:
	err = VecView(stateVarVecGlobal, PETSC_VIEWER_STDOUT_WORLD); //:DEBUG:
#endif

    // Copy current state vars and solution into superDM space
	PetscVec solutionVecGlobal = NULL, solutionVecLocal = solution.localVector();
	err = DMCreateGlobalVector(solutionDM, &solutionVecGlobal);
#if 1
	std::cout << std::endl << "solutionVecLocal" << std::endl; //:DEBUG:
	err = VecView(solutionVecLocal, PETSC_VIEWER_STDOUT_WORLD); //:DEBUG:
#endif
    err = DMLocalToGlobalBegin(solutionDM, solutionVecLocal, INSERT_VALUES, solutionVecGlobal);PYLITH_CHECK_ERROR(err);
    err = DMLocalToGlobalEnd(solutionDM, solutionVecLocal, INSERT_VALUES, solutionVecGlobal);PYLITH_CHECK_ERROR(err);
#if 1
	std::cout << std::endl << "solutionVecGlobal" << std::endl; //:DEBUG:
	err = VecView(solutionVecGlobal, PETSC_VIEWER_STDOUT_WORLD); //:DEBUG:
#endif
    err = VecISCopy(stateVarsSolutionGlobal, superIS[0], SCATTER_FORWARD, stateVarVecGlobal);PYLITH_CHECK_ERROR(err);
    err = VecISCopy(stateVarsSolutionGlobal, superIS[1], SCATTER_FORWARD, solutionVecGlobal);PYLITH_CHECK_ERROR(err);
#if 1
	PetscSection section = NULL;
	std::cout << std::endl << "stateVarsSolutionGlobal" << std::endl; //:DEBUG:
	err = VecView(stateVarsSolutionGlobal, PETSC_VIEWER_STDOUT_WORLD); //:DEBUG:
	std::cout << std::endl << "superDMSection" << std::endl; //:DEBUG:
	err = DMGetDefaultSection(superDM, &section); PYLITH_CHECK_ERROR(err); //:DEBUG:
	err = PetscSectionView(section, PETSC_VIEWER_STDOUT_WORLD); PYLITH_CHECK_ERROR(err); //:DEBUG:
#endif

	// Move superDM data to a local vector
    err = DMGlobalToLocalBegin(superDM, stateVarsSolutionGlobal, INSERT_VALUES, stateVarsSolutionLocal);PYLITH_CHECK_ERROR(err);
    err = DMGlobalToLocalEnd(superDM, stateVarsSolutionGlobal, INSERT_VALUES, stateVarsSolutionLocal);PYLITH_CHECK_ERROR(err);

    // Attach superDM data as auxiliary for updateStateVar kernel
    err = PetscObjectCompose((PetscObject) auxDM, "dmAux", (PetscObject) superDM);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) auxDM, "A",     (PetscObject) stateVarsSolutionLocal);PYLITH_CHECK_ERROR(err);
#else
    PetscDM dmState = _auxField->dmMesh();
    err = PetscObjectCompose((PetscObject) dmState, "dmAux", (PetscObject) solution.dmMesh());PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmState, "A", (PetscObject) solution.localVector());PYLITH_CHECK_ERROR(err);
#endif
    
    _setFEConstants(*_auxField, dt);

    // Set update kernel for each auxiliary subfield.
    PetscPointFunc* stateVarsKernels = (numSubfields > 0) ? new PetscPointFunc[numSubfields] : NULL;
    // By default, set all auxiliary subfield update kernels to NULL.
    for (size_t i = 0; i < numSubfields; ++i) {
        stateVarsKernels[i] = NULL;
    } // for
    for (UpdateStateVarsMap::const_iterator iter = _updateStateVarsKernels.begin(); iter != _updateStateVarsKernels.end(); ++iter) {
        const pylith::topology::Field::SubfieldInfo& sinfo = _auxField->subfieldInfo(iter->first.c_str());
        stateVarsKernels[sinfo.index] = iter->second;
    } // for

    err = DMProjectFieldLocal(auxDM, t, _auxField->localVector(), stateVarsKernels, INSERT_VALUES, _auxField->localVector());PYLITH_CHECK_ERROR(err);
    delete[] stateVarsKernels; stateVarsKernels = NULL;

#if 1
    // HAPPEN ONCE
    // Destroy superDM stuff
    for (size_t i = 0; i < 2; ++i) {err = ISDestroy(&superIS[i]);PYLITH_CHECK_ERROR(err);}
    err = DMDestroy(&superDM);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&stateVarIS);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&stateVarDM);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&stateVarsSolutionLocal);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&stateVarsSolutionGlobal);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&stateVarVecGlobal);PYLITH_CHECK_ERROR(err);
#endif
    
    PYLITH_METHOD_END;
} // updateStateVars


// ----------------------------------------------------------------------
// Write information (auxiliary field) output.
void
pylith::feassemble::IntegratorPointwise::writeInfo(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("writeInfo(void)");

    if (_output) {
        assert(_auxField);
        _output->writeInfo(*_auxField);
    } // if

    PYLITH_METHOD_END;
} // writeInfo


// ----------------------------------------------------------------------
// Write solution related output.
void
pylith::feassemble::IntegratorPointwise::writeTimeStep(const PylithReal t,
                                                       const PylithInt tindex,
                                                       const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("writeTimeStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<")");

    if (_output) {
        assert(_auxField);
        _output->writeTimeStep(t, tindex, solution, *_auxField);
    } // if

    PYLITH_METHOD_END;
} // writeTimeStep


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


// End of file
