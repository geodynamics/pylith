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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "DataWriter.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriter::DataWriter(void) :
  _timeScale(1.0),
  _numTimeSteps(0),
  _context("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriter::~DataWriter(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::DataWriter::deallocate(void)
{ // deallocate
} // deallocate
  

// ----------------------------------------------------------------------
// Set time scale for simulation time.
void
pylith::meshio::DataWriter::timeScale(const PylithScalar value)
{ // timeScale
  PYLITH_METHOD_BEGIN;

  if (value <= 0.0) {
    std::ostringstream msg;
    msg << "Time scale for simulation time (" << value << " must be positive.";
    throw std::runtime_error(msg.str());
  } // if
  
  _timeScale = value;

  PYLITH_METHOD_END;
} // timeScale
  
// ----------------------------------------------------------------------
// Prepare for writing files.
void
pylith::meshio::DataWriter::open(const topology::Mesh& mesh,
				 const int numTimeSteps,
				 const char* label,
				 const int labelId)
{ // open
  PYLITH_METHOD_BEGIN;

  _numTimeSteps = numTimeSteps;

  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  const char* meshName = NULL;
  PetscObjectGetName((PetscObject) dmMesh, &meshName);
  
  std::ostringstream s;
  s << "output_"
    << meshName;
  if (label)
    s << "_" << label << labelId;
  _context = s.str();

  PYLITH_METHOD_END;
} // open

// ----------------------------------------------------------------------
// Close output files.
void
pylith::meshio::DataWriter::close(void)
{ // close
  _context = "";
} // close

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
void
pylith::meshio::DataWriter::openTimeStep(const PylithScalar t,
					 const topology::Mesh& mesh,
					 const char* label,
					 const int labelId)
{ // openTimeStep
  // Default: no implementation.
} // openTimeStep

// ----------------------------------------------------------------------
// Cleanup after writing data for a time step.
void
pylith::meshio::DataWriter::closeTimeStep(void)
{ // closeTimeStep
  // Default: no implementation.
} // closeTimeStep

// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::DataWriter::DataWriter(const DataWriter& w) :
  _numTimeSteps(w._numTimeSteps),
  _context(w._context)
{ // copy constructor
} // copy constructor


// ----------------------------------------------------------------------
// Write dataset with names of points to file.
void
pylith::meshio::DataWriter::writePointNames(const char* const* names,
					    const int numNames)
{ // writePointNames
  // Default: no implementation.
} // writePointNames

// End of file 
