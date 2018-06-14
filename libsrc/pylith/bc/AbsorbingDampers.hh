// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/AbsorbingDampers.hh
 *
 * @brief C++ implementation of AbsorbingDampers for absorbing boundary conditions in dynamic problems.
 */

#if !defined(pylith_bc_absorbingdampers_hh)
#define pylith_bc_absorbingdampers_hh

// Include directives ---------------------------------------------------
#include "pylith/bc/bcfwd.hh" // forward delcaration
#include "pylith/feassemble/IntegratorBoundary.hh" // ISA IntegratorBoundary

#include "pylith/topology/topologyfwd.hh" // USES Field

// AbsorbingDampers ----------------------------------------------------
/// @brief AbsorbingDampers (e.g., traction) boundary conditions.
class pylith::bc::AbsorbingDampers : public pylith::feassemble::IntegratorBoundary {
    friend class TestAbsorbingDampers;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    AbsorbingDampers(void);

    /// Destructor.
    ~AbsorbingDampers(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * Create subfields in auxiliary fields (includes name of the field,
     * vector field type, discretization, and scale for
     * nondimensionalization) and set query functions for filling them
     * from a spatial database.
     *
     * @attention The order of the calls to subfieldAdd() must match the
     * order of the auxiliary fields in the FE kernels.
     *
     * @param[in] solution Solution field.
     */
    void _auxFieldSetup(const pylith::topology::Field& solution);

    /** Get factory for setting up auxliary fields.
     *
     * @returns Factor for auxiliary fields.
     */
    pylith::feassemble::AuxiliaryFactory* _auxFactory(void);

    /** Has point-wise functions (kernels) for integration/projection?
     *
     * @param[in] kernelsKey Set of kernels.
     * @returns True if we have kernels for that operation, otherwise false.
     */
    bool _hasFEKernels(const pylith::feassemble::IntegratorPointwise::FEKernelKeys kernelsKey) const;

    /** Set point-wise functions (kernels) for integration/projection.
     *
     * @param[in] solution Solution field.
     * @param[in] kernelsKey Set of kernels.
     */
    void _setFEKernels(const pylith::topology::Field& solution,
                       const pylith::feassemble::IntegratorPointwise::FEKernelKeys kernelsKey) const;


    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    pylith::topology::FieldBase::Description _description; ///< Description of field associated with BC.
    pylith::bc::AbsorbingDampersAuxiliaryFactory* _auxAbsorbingDampersFactory; ///< Factory for auxiliary subfields.
    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    AbsorbingDampers(const AbsorbingDampers&); ///< Not implemented.
    const AbsorbingDampers& operator=(const AbsorbingDampers&); ///< Not implemented.

}; // class AbsorbingDampers

#endif // pylith_bc_absorbingdampers_hh


// End of file
