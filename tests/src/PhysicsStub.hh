// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file libsrc/problems/PhysicsStub.hh
 *
 * @brief Minimal C++ implementation of Physics to allow testing of basic Physics functionality and use of Physics
 * objects in other tests.
 */

#if !defined(pylith_problems_physicsstub_hh)
#define pylith_problems_physicsstub_hh

#include "pylith/testing/testingfwd.hh" // forward declarations

#include "pylith/problems/Physics.hh" // ISA Physics

class pylith::problems::PhysicsStub : public pylith::problems::Physics {
    friend class TestPhysicsStub; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    PhysicsStub(void);

    /// Destructor
    ~PhysicsStub(void);

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
    std::vector<pylith::feassemble::Constraint*> createConstraints(const pylith::topology::Field& solution);

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ physicsMesh Finite-element mesh associated with physics.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& physicsMesh);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in\ physicsMesh Finite-element mesh associated with physics.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& physicsMesh);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::feassemble::AuxiliaryFactory* _auxiliaryFactory; ///< Factory for creating auxiliary subfields.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PhysicsStub(const PhysicsStub&); ///< Not implemented.
    const PhysicsStub& operator=(const PhysicsStub&); ///< Not implemented

}; // PhysicsStub

#endif // pylith_problems_physicsstub_hh

// End of file
