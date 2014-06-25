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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TimeDependent.hh" // implementation of object methods

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::TimeDependent::TimeDependent(void) :
  _dbInitial(0),
  _dbRate(0),
  _dbChange(0),
  _dbTimeHistory(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::TimeDependent::~TimeDependent(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::TimeDependent::deallocate(void)
{ // deallocate
  _dbInitial = 0; // TODO: Use shared pointers
  _dbRate = 0; // TODO: Use shared pointers
  _dbChange = 0; // TODO: Use shared pointers
  _dbTimeHistory = 0; // TODO: Use shared pointers
} // deallocate
  
// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::TimeDependent::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  if (!_dbChange && _dbTimeHistory) {
    std::ostringstream msg;
    msg << "Time dependent boundary condition '" << _getLabel() << "', has a "
	<< "time history database but not change in value spatial database.";
    throw std::runtime_error(msg.str());
  } // if
} // verifyConfiguration


// End of file 
