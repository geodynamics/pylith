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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesiveKin.hh" // implementation of object methods

#include "pylith/faults/KinSrc.hh" // USES KinSrc
#include "pylith/faults/AuxiliaryFactoryKinematic.hh" // USES AuxiliaryFactoryKinematic
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface
#include "pylith/feassemble/ConstraintSimple.hh" // USES ConstraintSimple

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::checkDiscretization()
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/fekernels/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensionalizer
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error
#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
// Return the number of impusles
void
pylith::faults::FaultCohesiveKinImpulses::getNumImpulses() {
  PetscErrorCode err = VecGetSize(_slipVecTotal, &N);PYLITH_CHECK_ERROR(err);
  return N;
} // getNumImpulses

// ---------------------------------------------------------------------------------------------------------------------
// Update slip subfield in auxiliary field at beginning of time step.
void
pylith::faults::FaultCohesiveKinImpulses::_updateSlip(pylith::topology::Field* auxiliaryField,
                                                      const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateSlip(auxiliaryField="<<auxiliaryField<<", t="<<t<<")");

    assert(auxiliaryField);
    assert(_normalizer);

    // Update slip subfield at current time step
    const PetscInt step = (PetscInt) t;
    PetscErrorCode err;

    if (step == 0) {
        // Create the amplitude field before the first impulse
        err = VecSet(_slipVecTotal, 0.0);PYLITH_CHECK_ERROR(err);
        const srcs_type::const_iterator rupturesEnd = _ruptures.end();
        for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter) {
            err = VecSet(_slipVecRupture, 0.0);PYLITH_CHECK_ERROR(err);

            KinSrc* src = r_iter->second;assert(src);
            src->updateSlip(_slipVecRupture, auxiliaryField, t, _normalizer->getTimeScale());
            err = VecAYPX(_slipVecTotal, 1.0, _slipVecRupture);
        } // for
    }

    // Transfer slip values from local PETSc slip vector to fault auxiliary field.
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryField->localSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

    pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField, "slip");
    PylithScalar* auxiliaryArray = auxiliaryVisitor.localArray();

    const PylithScalar* slipArray = NULL;
    err = VecGetArrayRead(_slipVecTotal, &slipArray);PYLITH_CHECK_ERROR(err);
    for (PetscInt p = pStart, iSlip = 0; p < pEnd; ++p) {
        const PetscInt slipDof = auxiliaryVisitor.sectionDof(p);
        const PetscInt slipOff = auxiliaryVisitor.sectionOffset(p);
        for (PetscInt iDof = 0; iDof < slipDof; ++iDof, ++iSlip) {
            auxiliaryArray[slipOff+iDof] = iSlip == step ? slipArray[iSlip] : 0.0;
        } // for
    } // for
    err = VecRestoreArrayRead(_slipVecTotal, &slipArray);PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        auxiliaryField->view("Fault auxiliary field after setting impulse.");
    } // if

    PYLITH_METHOD_END;
} // _updateSlip

// ---------------------------------------------------------------------------------------------------------------------

// End of file
