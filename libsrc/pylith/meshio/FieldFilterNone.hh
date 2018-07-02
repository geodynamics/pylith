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
 * @file libsrc/meshio/FieldFilterNone.hh
 *
 * @brief C++ object for no filtering.
 */

#if !defined(pylith_meshio_fieldfilternone_hh)
#define pylith_meshio_fieldfilternone_hh

// Include directives ---------------------------------------------------
#include "FieldFilter.hh" // ISA FieldFilter

// FieldFilterNone --------------------------------------------------
/** @brief C++ object for no filtering.
 */
class pylith::meshio::FieldFilterNone : public FieldFilter {

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    FieldFilterNone(void);

    /// Destructor
    ~FieldFilterNone(void);

    /** Create copy of filter.
     *
     * @returns Copy of filter.
     */
    FieldFilter* clone(void) const;

    /** Filter field.
     *
     * @param fieldIn Field to filter.
     * @returns Field field passed in without changes.
     */
    pylith::topology::Field* filter(pylith::topology::Field* fieldIn);

    // PROTECTED METHODS ////////////////////////////////////////////////////
protected:

    /** Copy constructor.
     *
     * @param f Filter to copy.
     * @returns Pointer to this.
     */
    FieldFilterNone(const FieldFilterNone& f);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    /// Not implemented.
    const FieldFilterNone& operator=(const FieldFilterNone&);

}; // FieldFilterNone

#endif // pylith_meshio_fieldfilternone_hh


// End of file
