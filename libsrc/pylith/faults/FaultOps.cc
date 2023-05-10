// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "FaultOps.hh" // implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory
#include "pylith/feassemble/IntegrationData.hh" // USES IntegrationData

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
// Create weighting vector for dynamic prescribed slip DAE.
void
pylith::faults::FaultOps::createDAEMassWeighting(pylith::feassemble::IntegrationData* integrationData) {
    PYLITH_METHOD_BEGIN;

    assert(integrationData);
    const pylith::topology::Field* solution = integrationData->getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const pylith::topology::Field::SubfieldInfo& velocityInfo = solution->getSubfieldInfo("velocity");
    const pylith::topology::Field::SubfieldInfo& lagrangeInfo = solution->getSubfieldInfo("lagrange_multiplier_fault");

    pylith::topology::Mesh* wtMesh = solution->getMesh().clone();
    pylith::topology::Field* wtField = new pylith::topology::Field(*wtMesh);
    wtField->setName("dae_mass_weighting");
    wtField->subfieldAdd(velocityInfo.description, velocityInfo.fe);
    wtField->subfieldAdd(lagrangeInfo.description, lagrangeInfo.fe);
    wtField->subfieldsSetup();
    wtField->createDiscretization();
    wtField->allocate();

    const char* dae_mass_weighting = pylith::feassemble::IntegrationData::dae_mass_weighting.c_str();
    integrationData->setMesh(dae_mass_weighting, wtMesh);
    integrationData->setField(dae_mass_weighting, wtField);

#if 0
    { // TEMPORARY DEBUGGING
        solution->view("SOLUTION", pylith::topology::Field::VIEW_LAYOUT);
        wtField->view("DAE WEIGHTING", pylith::topology::Field::VIEW_LAYOUT);
    } // TEMPORARY DEBUGGING
#endif

    PYLITH_METHOD_END;
} // createDAEMassWeighting


// ------------------------------------------------------------------------------------------------
// Create weighting vector for dynamic prescribed slip DAE.
void
pylith::faults::FaultOps::updateDAEMassWeighting(pylith::feassemble::IntegrationData* integrationData) {
    PYLITH_METHOD_BEGIN;
    PetscErrorCode err;

    const pylith::topology::Field* jacobianLumpedInv =
        integrationData->getField(pylith::feassemble::IntegrationData::lumped_jacobian_inverse);assert(jacobianLumpedInv);
    const pylith::topology::Field* weighting =
        integrationData->getField(pylith::feassemble::IntegrationData::dae_mass_weighting);assert(jacobianLumpedInv);

    const char* velocity = "velocity";
    pylith::topology::VecVisitorMesh domainVisitor(*jacobianLumpedInv, velocity);
    pylith::topology::VecVisitorMesh faultsVisitor(*weighting, velocity);
    PetscScalar* domainArray = domainVisitor.localArray();
    PetscScalar* faultsArray = faultsVisitor.localArray();

    err = VecSet(weighting->getLocalVector(), 1.0);PYLITH_CHECK_ERROR(err);
    PetscInt pStart = 0;
    PetscInt pEnd = 0;
    err = PetscSectionGetChart(faultsVisitor.selectedSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    for (PetscInt point = pStart; point < pEnd; ++point) {
        PetscInt numDof = faultsVisitor.sectionDof(point);
        if (numDof > 0) {
            assert(domainVisitor.sectionDof(point) == numDof);
            const PetscInt domainOff = domainVisitor.sectionOffset(point);
            const PetscInt faultsOff = faultsVisitor.sectionOffset(point);
            for (PetscInt iDof = 0; iDof < numDof; ++iDof) {
                faultsArray[faultsOff+iDof] = domainArray[domainOff+iDof];
            } // for

        } // if
    } // for

    { // TEMPORARY DEBUGGING
        weighting->view("DAE WEIGHTING", pylith::topology::Field::VIEW_ALL);
    } // TEMPORARY DEBUGGING

    PYLITH_METHOD_END;
} // updateDAEMassWeighting


// End of file
