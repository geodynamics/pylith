// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "GMVFile.hh" // implementation of class methods

#include "GMVFileAscii.hh"

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include <fstream> // uses std::fstream
#include <sstream> // uses std::ostringstream
#include <cstring> // uses strcmp()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------
pylith::meshio::GMVFile::GMVFile(const char* filename) :
  _filename(filename)
{ // constructor
} // constructor

// ----------------------------------------------------------------
pylith::meshio::GMVFile::~GMVFile(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------
bool
pylith::meshio::GMVFile::isAscii(const char* filename)
{ // isAscii
  PYLITH_METHOD_BEGIN;

  std::ifstream fin(filename);
  if (!(fin.is_open() && fin.good())) {
    std::ostringstream msg;
    msg << "Could not open GMV file '" << filename << "' for reading.";
    throw std::runtime_error(msg.str());
  } // if
  const int headerLen = strlen(GMVFileAscii::header())+1;
  char buffer[headerLen];
  fin.get(buffer, headerLen, '\n');
  fin.close();

  const bool result = (0 == strcmp(GMVFileAscii::header(), buffer)) ? true : false;
  PYLITH_METHOD_RETURN(result);
} // isAscii


// End of file 
