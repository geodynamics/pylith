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
 * @file libsrc/meshio/FieldFilter.hh
 *
 * @brief C++ object for filtering field upon output, etc.
 */

#if !defined(pylith_meshio_fieldfilter_hh)
#define pylith_meshio_fieldfilter_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field

// FieldFilter ---------------------------------------------------------
/** @brief C++ object for filtering field upon output, etc.
 */
class pylith::meshio::FieldFilter : public pylith::utils::PyreComponent {

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    FieldFilter(void);

    /// Destructor
    virtual ~FieldFilter(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Create copy of filter.
     *
     * @returns Copy of filter.
     */
    virtual
    FieldFilter* clone(void) const = 0;

    /** Filter field.
     *
     * @param fieldIn Field to filter.
     * @returns Field after applying filter.
     */
    virtual
    pylith::topology::Field* filter(pylith::topology::Field* fieldIn) = 0;

    // PROTECTED METHODS ////////////////////////////////////////////////////
protected:

    /** Copy constructor.
     *
     * @param f Filter to copy.
     * @returns Pointer to this.
     */
    FieldFilter(const FieldFilter& f);

private:
    /** operator=.
     *
     * @param f Filter to copy.
     * @returns Copy of filter.
     */
    const FieldFilter& operator=(const FieldFilter& f);

}; // FieldFilter

#endif // pylith_meshio_fieldfilter_hh


// End of file
