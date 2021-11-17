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

#include "InterfacePatches.hh" // implementation of object methods

#include "pylith/faults/TopologyOps.hh" // USES TopologyOps
#include "pylith/faults/FaultCohesive.hh" // USES FaultCohesive
#include "pylith/feassemble/FEKernelKey.hh" // USES FEKernelKey
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <utility> // USES std::pair
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

namespace pylith {
    namespace feassemble {
        namespace _InterfacePatches {
            inline
            PetscWeakForm
            getCellWeakForm(PetscDM dm,
                            const PetscInt cell) {
                PetscErrorCode err = 0;
                PetscDS prob = NULL;
                err = DMGetCellDS(dm, cell, &prob);PYLITH_CHECK_ERROR(err);
                PetscWeakForm weakForm = NULL;
                err = PetscDSGetWeakForm(prob, &weakForm);PYLITH_CHECK_ERROR(err);
                return weakForm;
            } // getCellWeakForm


        } // _InterfacePatches
    } // feassemble
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::InterfacePatches::InterfacePatches(void) :
    _labelName(pylith::topology::Mesh::getCellsLabelName()) {}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::InterfacePatches::~InterfacePatches(void) {}


// ------------------------------------------------------------------------------------------------
// Get name of label identifying integration patches.
const char*
pylith::feassemble::InterfacePatches::getLabelName(void) const {
    return _labelName.c_str();
} // getLabelName


// ------------------------------------------------------------------------------------------------
// Get weak form keys.
const pylith::feassemble::InterfacePatches::keysmap_t&
pylith::feassemble::InterfacePatches::getKeys(void) const {
    return _keys;
} // getKeys


// ------------------------------------------------------------------------------------------------
// Create integration patches corresponding to pairs of materials on negative and positive
pylith::feassemble::InterfacePatches*
pylith::feassemble::InterfacePatches::createMaterialPairs(const pylith::faults::FaultCohesive* fault,
                                                          const PetscDM dmSoln) {
    PYLITH_METHOD_BEGIN;

    assert(fault);
    const char* const cellsLabelName = pylith::topology::Mesh::getCellsLabelName();
    const std::string& patchLabelName = fault->getSurfaceMarkerLabel() + std::string("-integration-patches");
    InterfacePatches* patches = new InterfacePatches();assert(patches);
    patches->_labelName = patchLabelName;

    PetscErrorCode err = DMCreateLabel(dmSoln, patchLabelName.c_str());PYLITH_CHECK_ERROR(err);

    std::map<std::pair<int,int>, int> integrationPatches;
    PylithInt patchLabelValue = 0;

    PetscIS cohesiveCellsIS = NULL;
    PylithInt numCohesiveCells = 0;
    const PylithInt* cohesiveCells = NULL;
    err = DMGetStratumIS(dmSoln, cellsLabelName, fault->getInterfaceId(), &cohesiveCellsIS);PYLITH_CHECK_ERROR(err);
    err = ISGetSize(cohesiveCellsIS, &numCohesiveCells);PYLITH_CHECK_ERROR(err);assert(numCohesiveCells > 0);
    err = ISGetIndices(cohesiveCellsIS, &cohesiveCells);PYLITH_CHECK_ERROR(err);assert(cohesiveCells);

    for (PylithInt iCohesive = 0; iCohesive < numCohesiveCells; ++iCohesive) {
        const PetscInt cohesiveCell = cohesiveCells[iCohesive];
        assert(pylith::topology::MeshOps::isCohesiveCell(dmSoln, cohesiveCell));

        PetscInt adjacentCellNegative = -1;
        PetscInt adjacentCellPositive = -1;
        pylith::faults::TopologyOps::getAdjacentCells(&adjacentCellNegative, &adjacentCellPositive, dmSoln, cohesiveCell);
        assert(adjacentCellNegative >= 0);
        assert(adjacentCellPositive >= 0);

        std::pair<int, int> matPair;
        err = DMGetLabelValue(dmSoln, cellsLabelName, adjacentCellNegative, &matPair.first);
        err = DMGetLabelValue(dmSoln, cellsLabelName, adjacentCellPositive, &matPair.second);
        if (0 == integrationPatches.count(matPair)) {
            integrationPatches[matPair] = ++patchLabelValue;
            pythia::journal::debug_t debug("interfacepatches");
            // debug.activate();
            debug << pythia::journal::at(__HERE__)
                  << "Creating integration patch on fault '" << fault->getSurfaceMarkerLabel()
                  << "' for material pair ("<< matPair.first << "," << matPair.second<< ") "
                  << "using label '" << patchLabelName << "' with value " << patchLabelValue << "."
                  << pythia::journal::endl;

            // Get weak forms.
            PetscWeakForm weakFormCohesive = _InterfacePatches::getCellWeakForm(dmSoln, cohesiveCell);
            PetscWeakForm weakFormNegative = _InterfacePatches::getCellWeakForm(dmSoln, adjacentCellNegative);
            PetscWeakForm weakFormPositive = _InterfacePatches::getCellWeakForm(dmSoln, adjacentCellPositive);

            // Create keys
            WeakFormKeys weakFormKeys;
            pylith::feassemble::FEKernelKey* key = NULL;
            key = pylith::feassemble::FEKernelKey::create(weakFormCohesive, patchLabelName.c_str(), patchLabelValue);
            weakFormKeys.cohesive = *key;delete key;key = NULL;
            key = pylith::feassemble::FEKernelKey::create(weakFormNegative, cellsLabelName, matPair.first);
            weakFormKeys.negative = *key;delete key;key = NULL;
            key = pylith::feassemble::FEKernelKey::create(weakFormPositive, cellsLabelName, matPair.second);
            weakFormKeys.positive = *key;delete key;key = NULL;
            patches->_keys[patchLabelValue] = weakFormKeys;
        } // if
        err = DMSetLabelValue(dmSoln, patchLabelName.c_str(), cohesiveCell, integrationPatches[matPair]);
    } // for
    err = ISRestoreIndices(cohesiveCellsIS, &cohesiveCells);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&cohesiveCellsIS);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(patches);
} // createMaterialPairs


// End of file
