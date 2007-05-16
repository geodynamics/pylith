// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "Dirichlet.hh" // implementation of object methods

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::Dirichlet::Dirichlet(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::Dirichlet::~Dirichlet(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::Dirichlet::initialize(const ALE::Obj<ALE::Mesh>& mesh,
				  const spatialdata::geocoords::CoordSys* cs)
{ // initialize
} // initialize

// ----------------------------------------------------------------------
// Set constrained degrees of freedom in field.
void
pylith::bc::Dirichlet::setConstraints(const ALE::Obj<real_section_type>& field,
				      const ALE::Obj<ALE::Mesh>& mesh)
{ // setConstraints
} // setConstraints

// ----------------------------------------------------------------------
// Set values in field.
void
pylith::bc::Dirichlet::setField(const double t,
				const ALE::Obj<real_section_type>& field,
				const ALE::Obj<ALE::Mesh>& mesh)
{ // setField
} // setField


// End of file 
