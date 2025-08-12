// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/initializers/MeshInsertInterfaces.hh" // implementation of class methods

#include "pylith/problems/Problem.hh" // USES Problem
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Distributor.hh" // USES Distributor
#include "pylith/materials/Material.hh" // USES Material
#include "pylith/faults/FaultCohesive.hh" // USES FaultCohesive
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ------------------------------------------------------------------------------------------------
// Default constructor
pylith::initializers::MeshInsertInterfaces::MeshInsertInterfaces(void) {}


// ------------------------------------------------------------------------------------------------
// Default destructor
pylith::initializers::MeshInsertInterfaces::~MeshInsertInterfaces(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::initializers::MeshInsertInterfaces::deallocate(void) {
    InitializePhase::deallocate();
} // deallocate


#include <iostream>
// ------------------------------------------------------------------------------------------------
// Run initialization phase.
pylith::topology::Mesh*
pylith::initializers::MeshInsertInterfaces::run(pylith::topology::Mesh* mesh,
                                                const pylith::problems::Problem& problem) {
    PYLITH_METHOD_BEGIN;
    assert(mesh);

    if (!problem.getInterfaces().size()) {
        PetscErrorCode err = PETSC_SUCCESS;
        PetscDM dmOrig = mesh->getDM();assert(dmOrig);
        err = PetscObjectReference((PetscObject) dmOrig);PYLITH_CHECK_ERROR(err);
        pylith::topology::Mesh* meshNew = new pylith::topology::Mesh(dmOrig, *mesh);
        PYLITH_METHOD_RETURN(meshNew);
    } // if

    // Determine starting label value for cohesive cells.
    PylithInt cohesiveLabelValue = 100;
    for (auto material : problem.getMaterials()) {
        const PylithInt materialLabelValue = material->getLabelValue();
        cohesiveLabelValue = std::max(cohesiveLabelValue, materialLabelValue);
    } // for

    pylith::topology::Mesh* meshTmp = mesh;
    pylith::topology::Mesh* meshNew = nullptr;
    for (auto interface : problem.getInterfaces()) {
        interface->setCohesiveLabelValue(cohesiveLabelValue);
        meshNew = interface->transformTopology(meshTmp);
        if (meshTmp != mesh) {
            delete meshTmp;meshTmp = nullptr;
        } // if
        meshTmp = meshNew;
        cohesiveLabelValue += 1;
    } // for

    std::cout << "BEFORE OVERLAP" << std::endl;
    meshNew->view();
    meshNew->view(":before_overlap.tex:ascii_latex");

    PetscDM dmNew = nullptr;
    pylith::topology::Distributor::distributeOverlap(&dmNew, meshNew->getDM(), problem.getInterfaces());
    meshNew->setDM(dmNew);

    std::cout << std::endl << "AFTER OVERLAP" << std::endl;
    meshNew->view();
    meshNew->view(":after_overlap.tex:ascii_latex");

    PYLITH_METHOD_RETURN(meshNew);
} // run


// End of file
