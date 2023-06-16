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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/RheologyIncompressibleElasticity.hh
 *
 * @brief C++ abstract base class for bulk rheologies associated with the Incompressible elasticity equation.
 */

#if !defined(pylith_materials_rheologyIncompressibleelasticity_hh)
#define pylith_materials_rheologyIncompressibleelasticity_hh

#include "materialsfwd.hh"               // forward declarations
#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh"        // USES Field
#include "pylith/utils/arrayfwd.hh"              // USES std::vector
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain::ProjectKernels

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES Coordsys

#include "petscds.h" // USES PetscPointFunc, PetscPointJac

class pylith::materials::RheologyIncompressibleElasticity : public pylith::utils::PyreComponent
{
    friend class TestIsotropicLinearElasticity; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:
    /// Default constructor.
    RheologyIncompressibleElasticity(void);

    /// Destructor.
    virtual ~RheologyIncompressibleElasticity(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    virtual pylith::materials::AuxiliaryFactoryElasticity *getAuxiliaryFactory(void) = 0;

    /// Add rheology subfields to auxiliary field.
    virtual void addAuxiliarySubfields(void) = 0;

    /** Get f0p kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for potential.
     */
    virtual PetscPointFunc getKernelf0p(const spatialdata::geocoords::CoordSys *coordsys) const = 0;

    /** Get f0p kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for pressure.
     */
    virtual PetscPointFunc getKernelf1p(const spatialdata::geocoords::CoordSys *coordsys) const = 0;

    /** Get f1u kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual kernel for stress.
     */
    virtual PetscPointFunc getKernelf1u(const spatialdata::geocoords::CoordSys *coordsys) const = 0;

    /** Get Jf0pp kernel for LHS Jacobian F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jf0pp kernel.
     */
    virtual PetscPointJac getKernelJf0pp(const spatialdata::geocoords::CoordSys *coordsys) const = 0;

    /** Get Jf3uu kernel for LHS Jacobian F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jacobian kernel for elastic constants.
     */
    virtual PetscPointJac getKernelJf3uu(const spatialdata::geocoords::CoordSys *coordsys) const = 0;

    /** Get stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    virtual PetscPointFunc getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys *coordsys) const = 0;

    /** Add kernels for updating state variables.
     *
     * @param[inout] kernels Array of kernels for updating state variables.
     * @param[in] coordsys Coordinate system.
     */
    virtual void addKernelsUpdateStateVars(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels> *kernels,
                                           const spatialdata::geocoords::CoordSys *coordsys) const;

    /** Update kernel constants.
     *
     * @param[inout] kernelConstants Array of constants used in integration kernels.
     * @param[in] dt Current time step.
     */
    virtual void updateKernelConstants(pylith::real_array *kernelConstants,
                                       const PylithReal dt) const;

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:
    RheologyIncompressibleElasticity(const RheologyIncompressibleElasticity &);                  ///< Not implemented.
    const RheologyIncompressibleElasticity &operator=(const RheologyIncompressibleElasticity &); /// Not implemented.
};

// class RheologyIncompressibleElasticity

#endif // pylith_materials_rheologyIncompressibleelasticity_hh

// End of file
