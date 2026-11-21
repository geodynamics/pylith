// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


#pragma once

#include "pylith/feassemble/feassemblefwd.hh" // forward declarations

#include "pylith/feassemble/Integrator.hh" // USES EquationPart, eqnPart
#include "pylith/petsc/petsc_types.h" // HASA PetscDM
#include "pylith/topology/topologyfwd.hh" // USES Field

class pylith::feassemble::FEKernelKey {
    friend class TestFEKernelKey; // unit testing
    friend class TestInterfacePatches; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
  public:
    /// Default constructor.
    FEKernelKey(void);

    /// Default destructor.
    ~FEKernelKey(void);

    /** Factory for creating FEKernelKeyGet starting point.
     *
     * @param[in] weakForm PETSc weak form object.
     * @param[in] name Name of label designating integration domain.
     * @param[in] value Value of label designating integration domain.
     *
     * @return Key for finite-element integration.
     */
    static FEKernelKey* create(PetscWeakForm weakForm,
                               const char* name,
                               const int value);

    /** Get name of label.
     *
     * @returns Name of label.
     */
    const char* getName(void) const;

    /** Get value of label.
     *
     * @returns Label value.
     */
    int getValue(void) const;

    /** Get PETSc weak form.
     *
     * @returns PETSc weak form object.
     */
    PetscWeakForm getWeakForm(void) const;

    /** Get PETSc weak form key for integration.
     *
     * @param[in] solution Solution field.
     * @param[in] equationPart Equation part for weak form key.
     * @param[in] fieldTrial Name of solution subfield associated with trial function.
     * @param[in] fieldTrial Name of solution subfield associated with basis function.
     *
     * @returns PETSc weak form key.
     */
    PetscFormKey getPetscKey(const pylith::topology::Field& solution,
                             const PetscInt equationPart,
                             const char* fieldTrial = NULL,
                             const char* fieldBasis = NULL) const;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
  private:
    PetscWeakForm _weakForm; ///< PETSc weak form object associated with integration key.
    std::string _name; ///< Name of label designating integration domain.
    int _value; ///< Value of label designating integration domain.

}; // FEKernelKey

#include "FEKernelKey.icc"


// End of file
