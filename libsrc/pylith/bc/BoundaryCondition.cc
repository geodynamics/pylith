// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "BoundaryCondition.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh

#include <stdexcept> // USES std::runtime_error()

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::BoundaryCondition::BoundaryCondition(void) :
  _label("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::BoundaryCondition::~BoundaryCondition(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::BoundaryCondition::deallocate(void)
{ // deallocate
} // deallocate
  
// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::BoundaryCondition::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  throw std::logic_error("BoundaryCondition::verifyConfiguration(mesh) not implemented for PETSc dm.");
#if 0 // :MATT: Update this for PETSc dm.
  const ALE::Obj<topology::Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  if (!sieveMesh->hasIntSection(_label)) {
    std::ostringstream msg;
    msg << "Mesh missing group of vertices '" << _label
	<< "' for boundary condition.";
    throw std::runtime_error(msg.str());
  } // if
#endif
} // verifyConfiguration


// End of file 
