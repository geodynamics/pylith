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

/** @file libsrc/materials/Material.hh
 *
 * @brief C++ abstract base class for materials.
 */

#if !defined(pylith_materials_material_hh)
#define pylith_materials_material_hh

#include "materialsfwd.hh" // forward declarations

#include "pylith/problems/Physics.hh" // ISA Physics
#include "pylith/feassemble/IntegratorDomain.hh" // HOLDSA IntegratorDomain::ResidualKenels
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface::ResidualKernels

#include "pylith/utils/PetscOptions.hh" // USES PetscOptions

#include <string> // HASA std::string

// Material -------------------------------------------------------------
/** @brief C++ abstract base class for materials.
 *
 * Interface definition for a material. A material encapsulates both
 * the rheology as well as the governing equation.
 *
 * An individual material must abide by specific rules for the
 * interface, especially the order of the fields in the solution.
 *
 * Elasticity:
 *   + displacement, [velocity, Lagrange multipliers]
 *
 * Incompressible elasticity
 *   + displacement, pressure, [velocity, Lagrange multipliers]
 *
 * Poroelasticity
 *   + displacement, pressure, trace strain, [velocity, Lagrange multipliers]
 */

class pylith::materials::Material : public pylith::problems::Physics {
    friend class TestMaterial; // unit testing

    // PUBLIC TYPEDEFS ////////////////////////////////////////////////////////////////////////////
public:

    typedef pylith::feassemble::IntegratorInterface::ResidualKernels InterfaceResidualKernels;
    typedef pylith::feassemble::IntegratorInterface::JacobianKernels InterfaceJacobianKernels;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    Material(void);

    /// Destructor.
    virtual ~Material(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set descriptive label for material.
     *
     * @param value Label of material.
     */
    void setDescription(const char* value);

    /** Get descruptive label of material.
     *
     * @returns Label of material
     */
    const char* getDescription(void) const;

    /** Set gravity field.
     *
     * @param g Gravity field.
     */
    void setGravityField(spatialdata::spatialdb::GravityField* const g);

    /** Set body force kernels (intended for use in MMS tests).
     *
     * @param[in] kernels Array of body force kernels.
     */
    void setMMSBodyForceKernels(const std::vector<pylith::feassemble::IntegratorDomain::ResidualKernels>& kernels);

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    virtual
    std::vector<pylith::feassemble::Constraint*> createConstraints(const pylith::topology::Field& solution);

    /** Get default PETSc solver options appropriate for material.
     *
     * @param[in] isParallel True if running in parallel, False if running in serial.
     * @param[in] hasFault True if problem has fault, False otherwise.
     * @returns PETSc solver options.
     */
    virtual
    pylith::utils::PetscOptions* getSolverDefaults(const bool isParallel,
                                                   const bool hasFault) const;

    /** Get residual kernels for an interior interface bounding material.
     *
     * @param[in] solution Solution field.
     * @param[in] face Side of interior interface for kernels.
     * @returns Array of residual kernels for interior interface.
     */
    virtual
    std::vector<InterfaceResidualKernels> getInterfaceKernelsResidual(const pylith::topology::Field& solution,
                                                                      pylith::feassemble::IntegratorInterface::FaceEnum face) const;

    /** Get Jacobian kernels for an interior interface bounding material.
     *
     * @param[in] solution Solution field.
     * @param[in] face Side of interior interface for kernels.
     * @returns Array of Jacobian kernels for interior interface.
     */
    virtual
    std::vector<InterfaceJacobianKernels> getInterfaceKernelsJacobian(const pylith::topology::Field& solution,
                                                                      pylith::feassemble::IntegratorInterface::FaceEnum face) const;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field for gravitational body forces.
    std::vector<pylith::feassemble::IntegratorDomain::ResidualKernels> _mmsBodyForceKernels;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    std::string _description; ///< Descriptive label for material.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    Material(const Material&); ///< Not implemented.
    const Material& operator=(const Material&); ///< Not implemented

}; // Material

#endif // pylith_materials_material_hh

// End of file
