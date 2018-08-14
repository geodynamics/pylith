// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, Rice University
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2018 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/*** @file libsrc/feassemble/Physics.hh
 *
 * @brief Abstract base class for objects defining physics, such as behavior
 * of a bulk makterial, boundary condition, interface, or constraint.
 */

#if !defined(pylith_feassemble_physics_hh)
#define pylith_feassemble_physics_hh

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase
#include "pylith/utils/petscfwd.h" // USES PetscMat, PetscVec

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HASA spatialdb
#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional

class pylith::feassemble::Physics : public pylith::utils::PyreComponent {
    friend class TestPhysics; // unit testing

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /// Constructor
    Physics(void);

    /// Destructor
    virtual ~Physics(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Create integrator and set kernels.
     *
     * @solution Solution field.
     * @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

    /** Create constraint and set kernels.
     *
     * @solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution);

    /** Create auxiliary field.
     *
     * @mesh Finite-element mesh associated with physics.
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Mesh& mesh);

    /** Create derived field.
     *
     * @mesh Finite-element mesh associated with physics.
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createDerivedField(const pylith::topology::Mesh& mesh);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution) const = 0;

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalizer.

    pylith::utils::EventLogger* _logger; ///< Event logger.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    Physics(const Physics&); ///< Not implemented.
    const Physics& operator=(const Physics&); ///< Not implemented.

}; // Physics

#endif // pylith_feassemble_physics_hh

// End of file
