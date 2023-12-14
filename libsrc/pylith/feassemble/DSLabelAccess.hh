// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/feassemble/feassemblefwd.hh" // forward declarations

#include "pylith/utils/types.hh" // USES PetscInt
#include "pylith/utils/petscfwd.h" // HASA PetscDM, PetscDS, PetscWeakForm

#include <string> // HASA std::string

class pylith::feassemble::DSLabelAccess {
    friend class TestDSLabelAccess; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Constructor with PetscDM, label name, and label value.
     *
     * @param[in] dm PETSc DM holding label.
     * @param[in] labelName Name of label.
     * @param[in] labelValue Label value.
     */
    DSLabelAccess(const PetscDM dm,
                  const char* labelName,
                  const int labelValue);

    /// Default destructor
    ~DSLabelAccess(void);

    /** Get PETSc DM.
     *
     * @returns PETSc DM.
     */
    PetscDM dm(void) const;

    /** Get PETSc DMLabel.
     *
     * @returns PETSc DMLabel.
     */
    PetscDMLabel label(void) const;

    /** Get PETSc DMLabel value.
     *
     * @returns Label value.
     */
    PetscInt value(void) const;

    /** Get PetscDS.
     *
     * @returns PetscDS.
     */
    PetscDS ds(void) const;

    /** Get weak form associated with PetscDS.
     *
     * @returns Weak form.
     */
    PetscWeakForm weakForm(void) const;

    /** Get PETSc IS associated with label and value.
     *
     * @returns PETSc IS.
     */
    PetscIS cellsIS(void) const;

    /** Get number of cells in PETSc IS associated with label and value.
     *
     * @returns Number of cells.
     */
    PetscInt numCells(void) const;

    /// Remove overlap from list of cells in label.
    void removeOverlap(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    PetscDM _dm; ///< Cached PETSc dm for mesh.
    PetscDMLabel _label; ///< Cached PETSc DMLabel.
    PetscDS _ds; ///< Cached PetscDS for label and value.
    PetscWeakForm _weakForm; ///< Cached PETSc weak form for PetscDS.
    PetscIS _cellsIS; ///< Cached PETSc IS of cells for label and value.
    PetscInt _numCells; ///< Number of cells in PETSc IS.
    PetscInt _value; ///< Label value.
    std::string _name; ///< Name of label;

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    DSLabelAccess(void); ///< Not implemented.
    DSLabelAccess(const DSLabelAccess&); ///< Not implemented.
    const DSLabelAccess& operator=(const DSLabelAccess&); ///< Not implemented

};

#include "DSLabelAccess.icc"

// End of file
