// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "AbsorbingDampers.hh" // implementation of object methods

#include "AbsorbingDampersAuxiliaryFactory.hh" // USES AuxiliaryFactory

#include "pylith/fekernels/AbsorbingDampers.hh" // USES AbsorbingDampers kernels
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> \
    // USES std::ostringstream

// ----------------------------------------------------------------------
const char* pylith::bc::AbsorbingDampers::_pyreComponent = "absorbingdampers";

// Local "private" functions.
namespace pylith {
    namespace bc {
        static void _setFEKernelsRHSResidual(const AbsorbingDampers* const bc,
                                             PetscDS prob,
                                             const PylithInt fieldIndex);
    } // bc
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::AbsorbingDampers::AbsorbingDampers(void) :
    _auxAbsorbingDampersFactory(new pylith::bc::AbsorbingDampersAuxiliaryFactory)
{ // constructor
    PyreComponent::name(_pyreComponent);

    _field = "velocity";
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::AbsorbingDampers::~AbsorbingDampers(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::AbsorbingDampers::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    IntegratorBoundary::deallocate();

    delete _auxAbsorbingDampersFactory; _auxAbsorbingDampersFactory = NULL;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::AbsorbingDampers::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    IntegratorBoundary::verifyConfiguration(solution);

    const pylith::topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    if (pylith::topology::Field::VECTOR != info.description.vectorFieldType) {
        std::ostringstream msg;
        msg << "Absorbing boundary condition cannot be applied to non-vector field '"<< _field
            << "' in solution.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Setup auxiliary subfields (discretization and query fns).
void
pylith::bc::AbsorbingDampers::_auxFieldSetup(const pylith::topology::Field& solution)
{ // _auxFieldSetup
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxFieldSetup(solution="<<solution.label()<<")");

    assert(_auxAbsorbingDampersFactory);
    assert(_normalizer);
    _auxAbsorbingDampersFactory->initialize(_auxField, *_normalizer, solution.spaceDim(),
                                            &solution.subfieldInfo(_field.c_str()).description);

    // :ATTENTION: The order of the factory methods must match the order of the auxiliary subfields in the FE kernels.

    _auxAbsorbingDampersFactory->density();
    _auxAbsorbingDampersFactory->vp();
    _auxAbsorbingDampersFactory->vs();

    PYLITH_METHOD_END;
} // _auxFieldSetup


// ----------------------------------------------------------------------
// Get factory for setting up auxliary fields.
pylith::feassemble::AuxiliaryFactory*
pylith::bc::AbsorbingDampers::_auxFactory(void) {
    return _auxAbsorbingDampersFactory;
} // _auxFactory


// ----------------------------------------------------------------------
// Does boundary conditon have point-wise functions (kernels) for integration/projection.
bool
pylith::bc::AbsorbingDampers::_hasFEKernels(IntegratorPointwise::FEKernelKeys kernelsKey) const {
    bool hasKernels = false;
    switch (kernelsKey) {
    case KERNELS_RHS_RESIDUAL:
        hasKernels = true;
        break;
    case KERNELS_LHS_RESIDUAL:
    case KERNELS_RHS_JACOBIAN:
    case KERNELS_LHS_JACOBIAN:
    case KERNELS_LHS_JACOBIAN_LUMPEDINV:
    case KERNELS_UPDATE_STATE_VARS:
    case KERNELS_DERIVED_FIELDS:
        hasKernels = false;
        break;
    default:
        PYLITH_COMPONENT_ERROR("Unrecognized finite-element kernels key '"<<kernelsKey<<"'.");
        throw std::logic_error("Unrecognized finite-element kernels key.");
    } // switch
    return hasKernels;
} // _hasKernels


// ----------------------------------------------------------------------
// Set point-wise functions (kernels) for integration/projection.
void
pylith::bc::AbsorbingDampers::_setFEKernels(const pylith::topology::Field& solution,
                                            IntegratorPointwise::FEKernelKeys kernelsKey) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernels(solution="<<solution.label()<<", kernelsKey="<<kernelsKey<<")");

    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    const int fieldIndex = solution.subfieldInfo(_field.c_str()).index;

    switch (kernelsKey) {
    case KERNELS_RHS_RESIDUAL:
        _setFEKernelsRHSResidual(this, prob, fieldIndex);
        break;
    case KERNELS_LHS_RESIDUAL:
    case KERNELS_RHS_JACOBIAN:
    case KERNELS_LHS_JACOBIAN:
    case KERNELS_LHS_JACOBIAN_LUMPEDINV:
    case KERNELS_UPDATE_STATE_VARS:
    case KERNELS_DERIVED_FIELDS:
        break;
    } // switch
} // _setFEKernels


// ----------------------------------------------------------------------
// Set point-wise functions (kernels) for integrating RHS residual for vector field.
void
pylith::bc::_setFEKernelsRHSResidual(const AbsorbingDampers* const bc,
                                     PetscDS prob,
                                     const PylithInt fieldIndex) {
    PYLITH_METHOD_BEGIN;

    PetscBdPointFunc g0 = pylith::fekernels::AbsorbingDampers::g0;
    PetscBdPointFunc g1 = NULL;

    PetscErrorCode err = PetscDSSetBdResidual(prob, fieldIndex, g0, g1); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSResidualVector


// End of file
