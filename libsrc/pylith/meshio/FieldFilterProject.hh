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
 * @file libsrc/meshio/FieldFilterProject.hh
 *
 * @brief C++ object for projecting a field to another basis.
 */

#if !defined(pylith_meshio_fieldfilterproject_hh)
#define pylith_meshio_fieldfilterproject_hh

#include "FieldFilter.hh" // ISA FieldFilter

class pylith::meshio::FieldFilterProject : public FieldFilter {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    FieldFilterProject(void);

    /// Destructor
    ~FieldFilterProject(void);

    /** Create copy of filter.
     *
     * @returns Copy of filter.
     */
    FieldFilter* clone(void) const;

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set basis order for projected field.
     *
     * @param[in] value Basis order.
     */
    void basisOrder(const int value);

    /** Filter vertex field.
     *
     * @returns Field after applying filter.
     * @param fieldIn Field to filter.
     */
    pylith::topology::Field* filter(pylith::topology::Field* fieldIn);

    /** Pass thru solution in pointwise function.
     *
     * We pass through the solution to the resulting field. The auxiliary field is ignored.
     *
     * @param[in] dim Spatial dimension.
     * @param[in] numS Number of registered subfields in solution field.
     * @param[in] numA Number of registered subfields in auxiliary field.
     * @param[in] sOff Offset of registered subfields in solution field [numS].
     * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
     * @param[in] s Solution field with all subfields.
     * @param[in] s_t Time derivative of solution field.
     * @param[in] s_x Gradient of solution field.
     * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
     * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
     * @param[in] a Auxiliary field with all subfields.
     * @param[in] a_t Time derivative of auxiliary field.
     * @param[in] a_x Gradient of auxiliary field.
     * @param[in] t Time for residual evaluation.
     * @param[in] x Coordinates of point evaluation.
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] field [dim].
     */
    static
    void passThruSoln(const PylithInt dim,
                      const PylithInt numS,
                      const PylithInt numA,
                      const PylithInt sOff[],
                      const PylithInt sOff_x[],
                      const PylithScalar s[],
                      const PylithScalar s_t[],
                      const PylithScalar s_x[],
                      const PylithInt aOff[],
                      const PylithInt aOff_x[],
                      const PylithScalar a[],
                      const PylithScalar a_t[],
                      const PylithScalar a_x[],
                      const PylithReal t,
                      const PylithScalar x[],
                      const PylithInt numConstants,
                      const PylithScalar constants[],
                      PylithScalar field[]);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Copy constructor.
     *
     * @param f Filter to copy.
     * @returns Pointer to this.
     */
    FieldFilterProject(const FieldFilterProject& f);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::topology::Field* _fieldProj; ///< Projected field.
    PetscPointFunc* _passThruFns; ///< Pass through point functions.
    int _basisOrder; ///< Basis order for projected field.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /// Not implemented.
    const FieldFilterProject& operator=(const FieldFilterProject&);

}; // FieldFilterProject

#endif // pylith_meshio_fieldfilterproject_hh

// End of file
