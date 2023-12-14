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

/*
 * In some cases the weak form involves both integration over the cohesive cell as well as
 * the adjacent faces on the negative and positive sides of the fault. In such cases the pointwise
 * functions depend on the materials, so we split the integration into patches. Each patch
 * corresponds to a pair of materials on the negative and positive sides of the fault.
 */

#include "pylith/feassemble/feassemblefwd.hh" // forward declarations

#include "pylith/faults/faultsfwd.hh" // USES FaultCohesive

#include "pylith/feassemble/FEKernelKey.hh" // HASA FEKernelKey

#include <string> // HASA std::string
#include <map> // HASA std::map

class pylith::feassemble::InterfacePatches {
    friend class TestInterfacePatches; // unit testing

    // PUBLIC STRUCTS /////////////////////////////////////////////////////////////////////////////
public:

    struct WeakFormKeys {
        FEKernelKey cohesive;
        FEKernelKey negative;
        FEKernelKey positive;
    }; // WeakFormKeys

    typedef std::map<int, WeakFormKeys> keysmap_t;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    InterfacePatches(void);

    /// Destructor.
    virtual ~InterfacePatches(void);

    /** Get name of label identifying integration patches.
     *
     * @returns Name of label.
     */
    const char* getLabelName(void) const;

    /** Get weak form keys.
     *
     * @returns Weak form keys
     */
    const keysmap_t& getKeys(void) const;

    /** Create integration patches corresponding to pairs of materials on negative and positive
     * sides of the fault.
     *
     * @param[in] fault Fault associated with interface.
     * @param[in] dmSoln PETSc DM associated with solution.
     */
    static
    InterfacePatches* createMaterialPairs(const pylith::faults::FaultCohesive* fault,
                                          const PetscDM dmSoln);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    std::string _labelName; ///< Name of label identifying integration patches.
    keysmap_t _keys; ///< Weak form keys.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    InterfacePatches(const InterfacePatches&); ///< Not implemented
    const InterfacePatches& operator=(const InterfacePatches&); ///< Not implemented

}; // class InterfacePatches

// End of file
