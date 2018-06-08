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
 * @file libsrc/meshio/VertexFilterDecimateP1.hh
 *
 * @brief C++ object for decimating a vertex field to P1.
 */

#if !defined(pylith_meshio_vertexfilterdecimatep1_hh)
#define pylith_meshio_vertexfilterdecimatep1_hh

// Include directives ---------------------------------------------------
#include "VertexFilter.hh" // ISA VertexFilter

// VertexFilterDecimateP1 --------------------------------------------------
/** @brief C++ object for computing vector norms for fields over
 * vertices when outputing finite-element data.
 */
class pylith::meshio::VertexFilterDecimateP1 : public VertexFilter { // VertexFilterDecimateP1

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    VertexFilterDecimateP1(void);

    /// Destructor
    ~VertexFilterDecimateP1(void);

    /** Create copy of filter.
     *
     * @returns Copy of filter.
     */
    VertexFilter* clone(void) const;

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Filter vertex field.
     *
     * @param fieldIn Field to filter.
     */
    topology::Field& filter(const topology::Field& fieldIn);

    // PROTECTED METHODS ////////////////////////////////////////////////////
protected:

    /** Copy constructor.
     *
     * @param f Filter to copy.
     * @returns Pointer to this.
     */
    VertexFilterDecimateP1(const VertexFilterDecimateP1& f);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    /// Not implemented.
    const VertexFilterDecimateP1& operator=(const VertexFilterDecimateP1&);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    topology::Field* _fieldP1; ///< Filtered vertex field

}; // VertexFilterDecimateP1

#endif // pylith_meshio_vertexfilterdecimatep1_hh


// End of file
