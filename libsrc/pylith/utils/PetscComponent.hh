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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/utils/PetscComponent.hh
 *
 * @brief C++ abstract base class for managing objects that hold PETSc objects.
 */

#if !defined(pylith_utils_petsccomponent_hh)
#define pylith_utils_petsccomponent_hh

// Include directives ---------------------------------------------------
#include "utilsfwd.hh" // forward declarations

// PetscComponent ----------------------------------------------------------
/** @brief C++ abstract base class for managing objects that hold PETSc objects.
 *
 * Facilitates deallocation of PETSc objects before PetscFinalize().
 */
class pylith::utils::PetscComponent {
    friend class TestPetscComponent; // unit testing

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /// Constructor
    PetscComponent(void);

    /// Destructor
    virtual ~PetscComponent(void);

    /// Deallocate PETSc objects.
    virtual
    void deallocate(void) = 0;

}; // PetscComponent

#endif // pylith_utils_petsccomponent_hh


// End of file
