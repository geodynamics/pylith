// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/DSLabelAccess.hh
 *
 * @brief C++ helper class for accessing PetscDMLabel, PetscDS, and PetscWeakForm
 * for a given DM, label, and label value.
 */

#if !defined(pylith_feassemble_dslabelaccess_hh)
#define pylith_feassemble_dslabelaccess_hh

#include "feassemblefwd.hh" // forward declarations

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

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    PetscDM _dm; ///< Cached PETSc dm for mesh.
    PetscDMLabel _label; ///< Cached PETSc DMLabel.
    PetscDS _ds; ///< Cached PetscDS for label and value.
    PetscWeakForm _weakForm; ///< Cached PETSc weak form for PetscDS.
    PetscIS _cellsIS; ///< Cached PETSc IS of cells for label and value.
    PetscInt _numCells; ///< Number of cells in PETSc IS.
    PetscInt _value; ///< Label value.
    std::string _name; ///< Name of label;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    DSLabelAccess(const DSLabelAccess&); ///< Not implemented
    const DSLabelAccess& operator=(const DSLabelAccess&); ///< Not implemented

};

#include "DSLabelAccess.icc"

#endif // pylith_feassemble_dslabelaccess_hh

// End of file
