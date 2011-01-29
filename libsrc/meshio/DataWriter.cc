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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriter<mesh_type, field_type>::DataWriter(void) :
  _numTimeSteps(0),
  _context("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriter<mesh_type, field_type>::~DataWriter(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriter<mesh_type, field_type>::deallocate(void)
{ // deallocate
} // deallocate
  
// ----------------------------------------------------------------------
// Prepare for writing files.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriter<mesh_type, field_type>::open(const mesh_type& mesh,
					    const int numTimeSteps,
					    const char* label,
					    const int labelId)
{ // open
  _numTimeSteps = numTimeSteps;

  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  _context = std::string("output_") + sieveMesh->getName() + 
    ((label) ? std::string("_") + std::string(label) : std::string(""));
} // open

// ----------------------------------------------------------------------
// Close output files.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriter<mesh_type, field_type>::close(void)
{ // close
  _context = "";
} // close

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriter<mesh_type, field_type>::openTimeStep(const double t,
						    const mesh_type& mesh,
						    const char* label,
						    const int labelId)
{ // openTimeStep
} // openTimeStep

// ----------------------------------------------------------------------
// Cleanup after writing data for a time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriter<mesh_type, field_type>::closeTimeStep(void)
{ // closeTimeStep
} // closeTimeStep

// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriter<mesh_type, field_type>::DataWriter(const DataWriter& w) :
  _numTimeSteps(w._numTimeSteps),
  _context(w._context)
{ // copy constructor
} // copy constructor


// End of file 
