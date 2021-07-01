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

/** @file libsrc/feassemble/PhysicsImplementation.hh
 *
 * @brief C++ abstract base class defining interface for constraining degrees of freedom in the solution.
 */

#if !defined(pylith_feassemble_physicaimplementation_hh)
#define pylith_feassemble_physicaimplementation_hh

#include "pylith/feassemble/feassemblefwd.hh"

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/problems/problemsfwd.hh" // HASA Physics, ObserversPhysics
#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/utils/array.hh" // HASA int_array
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

class pylith::feassemble::PhysicsImplementation : public pylith::utils::GenericComponent {
    friend class TestPhysicsImplementation; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] physics Physics implemented by constraint.
     */
    PhysicsImplementation(pylith::problems::Physics* const physics);

    /// Destructor.
    virtual ~PhysicsImplementation(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Get mesh associated with domain governed by specified physics.
     *
     * @returns Mesh associated with constrained domain.
     */
    virtual
    const pylith::topology::Mesh& getPhysicsDomainMesh(void) const = 0;

    /** Get auxiliary field.
     *
     * @returns field Field over boundary.
     */
    const pylith::topology::Field* getAuxiliaryField(void) const;

    /** Get derived field.
     *
     * @return field Field over integrator domain.
     */
    const pylith::topology::Field* getDerivedField(void) const;

    /** Notify observers of current solution.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     */
    void notifyObservers(const PylithReal t,
                         const PylithInt tindex,
                         const pylith::topology::Field& solution);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::problems::Physics* const _physics; ///< Physics associated with constraint.
    pylith::topology::Field* _auxiliaryField; ///< Auxiliary field for this constraint.
    pylith::topology::Field* _derivedField; ///< Derived field for this constraint.
    pylith::problems::ObserversPhysics* _observers; ///< Observers component.

    pylith::utils::EventLogger* _logger; ///< Event logger.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PhysicsImplementation(void); /// Not implemented.
    PhysicsImplementation(const PhysicsImplementation &); ///< Not implemented
    const PhysicsImplementation& operator=(const PhysicsImplementation&); ///< Not implemented

}; // class PhysicsImplementation

#endif // pylith_feassemble_physicaimplementation_hh

// End of file
