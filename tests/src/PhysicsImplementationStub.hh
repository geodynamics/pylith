// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

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

// End of file
