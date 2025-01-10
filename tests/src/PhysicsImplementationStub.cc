// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "PhysicsImplementationStub.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "tests/src/StubMethodTracker.hh" // USES StubMethodTracker

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::feassemble::PhysicsImplementationStub::PhysicsImplementationStub(void) :
    PhysicsImplementation(NULL),
    _mesh(new pylith::topology::Mesh) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::feassemble::PhysicsImplementationStub::~PhysicsImplementationStub(void) {
    delete _mesh;_mesh = NULL;
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Get mesh associated with domain governed by specified physics.
const pylith::topology::Mesh&
pylith::feassemble::PhysicsImplementationStub::getPhysicsDomainMesh(void) const {
    pylith::testing::StubMethodTracker tracker("pylith::feassemble::PhysicsImplementationStub::getPhysicsDomainMesh");

    assert(_mesh);
    return *_mesh;
} // verifyConfiguration


// End of file
