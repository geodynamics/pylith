// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2019 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "PhysicsImplementationStub.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/testing/StubMethodTracker.hh" // USES StubMethodTracker

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
