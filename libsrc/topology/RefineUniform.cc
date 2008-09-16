// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "RefineUniform.hh" // implementation of class methods

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::topology::RefineUniform::RefineUniform(void)
{ // constructor
} // constructor
 
// ----------------------------------------------------------------------
// Destructor
pylith::topology::RefineUniform::~RefineUniform(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Refine mesh.
void
pylith::topology::RefineUniform::refine(ALE::Obj<Mesh>* const newMesh,
					const ALE::Obj<Mesh>& mesh,
					const int levels)
{ // refine
} // refine


// End of file 
