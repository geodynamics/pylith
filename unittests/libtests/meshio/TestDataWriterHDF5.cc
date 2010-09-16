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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterHDF5.hh" // Implementation of class methods

#include <string.h> // USES strcmp()
#include <iostream> // USES std::cerr
#include <fstream> // USES std::ifstream

// ----------------------------------------------------------------------
// Check HDF5 file against archived file.
void
pylith::meshio::TestDataWriterHDF5::checkFile(const char* filename)
{ // checkFile

  const std::string filenameE = "data/" + std::string(filename);

} // checkFile


// End of file 
