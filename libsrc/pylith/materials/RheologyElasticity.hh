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
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/RheologyElasticity.hh
 *
 * @brief C++ abstract base class for bulk rheologies associated with the elasticity equation.
 */

#if !defined(pylith_materials_rheologyelasticity_hh)
#define pylith_materials_rheologyelasticity_hh

#include "materialsfwd.hh" // forward declarations
#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/arrayfwd.hh" // USES std::vector
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain::ProjectKernels

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES Coordsys

#include "petscds.h" // USES PetscPointFunc, PetscPointJac

class pylith::materials::RheologyElasticity : public pylith::utils::PyreComponent {
    friend class TestIsotropicLinearElasticity; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    RheologyElasticity(void);

    /// Destructor.
    virtual ~RheologyElasticity(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    virtual
    pylith::materials::AuxiliaryFactoryElasticity* getAuxiliaryFactory(void) = 0;

    /// Add rheology subfields to auxiliary field.
    virtual
    void addAuxiliarySubfields(void) = 0;

    /** Get stress kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for stress.
     */
    virtual
    PetscPointFunc getKernelResidualStress(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get elastic constants kernel for RHS Jacobian G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS Jacobian kernel for elastic constants.
     */
    virtual
    PetscPointJac getKernelJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get f0 kernel for LHS interface residual, F(t,s,dot{s}), for negative fault face.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual f0 kernel.
     */
    virtual
    PetscBdPointFunc getInterfaceKernelResidualF0Neg(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get f0 kernel for LHS interface residual, F(t,s,dot{s}), for positive fault face.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual f0 kernel.
     */
    virtual
    PetscBdPointFunc getInterfaceKernelResidualF0Pos(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get f1 kernel for LHS interface residual, F(t,s,dot{s}), for negative fault face.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual f1 kernel.
     */
    virtual
    PetscBdPointFunc getInterfaceKernelResidualF1Neg(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get f1 kernel for LHS interface residual, F(t,s,dot{s}), for positive fault face.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual f1 kernel.
     */
    virtual
    PetscBdPointFunc getInterfaceKernelResidualF1Pos(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get Jf1lu kernel for LHS Jacobian F(t,s,dot{s}) for negative fault face.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jacobian kernel Jf1lu.
     */
    virtual
    PetscBdPointJac getInterfaceKernelJacobianF1Neg(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get Jf1lu kernel for LHS Jacobian F(t,s,dot{s}) for positive fault face.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jacobian kernel Jf1lu.
     */
    virtual
    PetscBdPointJac getInterfaceKernelJacobianF1Pos(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get Jf3lu kernel for LHS Jacobian F(t,s,dot{s}) for negative fault face.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jacobian kernel Jf3lu.
     */
    virtual
    PetscBdPointJac getInterfaceKernelJacobianF3Neg(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get Jf3lu kernel for LHS Jacobian F(t,s,dot{s}) for positive fault face.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jacobian kernel Jf3lu.
     */
    virtual
    PetscBdPointJac getInterfaceKernelJacobianF3Pos(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get triggers for needing to compute the elastic constants for the RHS Jacobian.
     *
     * @returns Triggers for needing to recompute the RHS Jacobian.
     */
    int getLHSJacobianTriggers(void) const;

    /** Get stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    virtual
    PetscPointFunc getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Add kernels for updating state variables.
     *
     * @param[inout] kernels Array of kernels for updating state variables.
     * @param[in] coordsys Coordinate system.
     */
    virtual
    void addKernelsUpdateStateVars(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                   const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Update kernel constants.
     *
     * @param[inout] kernelConstants Array of constants used in integration kernels.
     * @param[in] dt Current time step.
     */
    virtual
    void updateKernelConstants(pylith::real_array* kernelConstants,
                               const PylithReal dt) const;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////

    int _lhsJacobianTriggers; ///< Triggers for needing to recompute the RHS Jacobian.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    RheologyElasticity(const RheologyElasticity&); ///< Not implemented.
    const RheologyElasticity& operator=(const RheologyElasticity&); /// Not implemented.

}; // class RheologyElasticity

#endif // pylith_materials_rheologyelasticity_hh

// End of file
