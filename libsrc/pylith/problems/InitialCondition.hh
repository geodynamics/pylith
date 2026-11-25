// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


#pragma once

#include "pylith/problems/problemsfwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/utils/arrayfwd.hh" // USES string_vector

#include "pylith/scales/scalesfwd.hh" // HASA Scales
#include "pylith/topology/topologyfwd.hh" // USES Field

class pylith::problems::InitialCondition : public pylith::utils::PyreComponent {
    friend class TestInitialCondition; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
  public:
    /// Constructor
    InitialCondition(void);

    /// Destructor
    virtual ~InitialCondition(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set solution subfields for initial condition.
     *
     * @param[in] subfields Array of names of solution subfields.
     * @param[in] numSubfields Number of subfields.
     */
    void setSubfields(const char* subfields[],
                      const int numSubfields);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    virtual void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Set solution to values for initial condition.
     *
     * @param[out] solution Solution field.
     * @param[in] scales Nondimensionalization.
     */
    virtual void setValues(pylith::topology::Field* solution,
                           const pylith::scales::Scales& scales) = 0;

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
  protected:
    pylith::string_vector _subfields; ///< Names of fields for initial conditions.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
  private:
    InitialCondition(const InitialCondition&); ///< Not implemented
    const InitialCondition& operator=(const InitialCondition&); ///< Not implemented

}; // InitialCondition


// End of file
