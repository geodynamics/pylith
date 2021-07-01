// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/AbsorbingDampers.hh
 *
 * @brief C++ implementation of AbsorbingDampers for absorbing boundary conditions in dynamic problems.
 */

#if !defined(pylith_bc_absorbingdampers_hh)
#define pylith_bc_absorbingdampers_hh

#include "pylith/bc/BoundaryCondition.hh" // ISA Physics

#include "pylith/topology/topologyfwd.hh" // USES Field

class pylith::bc::AbsorbingDampers : public pylith::bc::BoundaryCondition {
    friend class TestAbsorbingDampers; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
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

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution);

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& domainMesh);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    void _updateKernelConstants(const PylithReal dt);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::bc::AbsorbingDampersAuxiliaryFactory* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    PylithReal _refDir1[3]; ///< First choice reference direction used to compute boundary tangential directions.
    PylithReal _refDir2[3]; ///< Second choice reference direction used to compute boundary tangential directions.
    std::string _boundaryLabel; ///< Label to identify boundary condition points in mesh.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AbsorbingDampers(const AbsorbingDampers&); ///< Not implemented.
    const AbsorbingDampers& operator=(const AbsorbingDampers&); ///< Not implemented.

};

// class AbsorbingDampers

#endif // pylith_bc_absorbingdampers_hh

// End of file
