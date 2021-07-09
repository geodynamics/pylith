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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/PhysicsImplementationStub.hh
 *
 * @brief Minimal C++ implementation of PhysicsImplementation to allow testing of basic PhysicsImplementation
 * functionality and use of PhysicsImplementation objects in other tests.
 */

#if !defined(pylith_feassemble_physicsimplementationstub_hh)
#define pylith_feassemble_physicsimplementationstub_hh

#include "pylith/testing/testingfwd.hh" // forward declarations

#include "pylith/feassemble/PhysicsImplementation.hh" // ISA PhysicsImplementation

#include "pylith/topology/topologyfwd.hh" // HOLDSA Mesh

class pylith::feassemble::PhysicsImplementationStub : public pylith::feassemble::PhysicsImplementation {
    friend class TestPhysicsImplementationStub; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    PhysicsImplementationStub(void);

    /// Destructor
    ~PhysicsImplementationStub(void);

    /** Get mesh associated with domain governed by specified physics.
     *
     * @returns Mesh associated with constrained domain.
     */
    const pylith::topology::Mesh& getPhysicsDomainMesh(void) const;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PhysicsImplementationStub(const PhysicsImplementationStub&); ///< Not implemented.
    const PhysicsImplementationStub& operator=(const PhysicsImplementationStub&); ///< Not implemented

    pylith::topology::Mesh* _mesh; ///< Mesh object.

}; // PhysicsImplementationStub

#endif // pylith_feassemble_physicsimplementationstub_hh

// End of file
