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

#include "SolutionIO.hh" // implementation of class methods

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::SolutionIO::SolutionIO(void) :
  _cs(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::SolutionIO::~SolutionIO(void)
{ // destructor
  delete _cs; _cs = 0;
} // destructor  

// ----------------------------------------------------------------------
// Set coordinate system for output.
void
pylith::meshio::SolutionIO::coordsys(const spatialdata::geocoords::CoordSys* cs)
{ // coordsys
  delete _cs; _cs = (0 != cs) ? cs->clone() : 0;
} // coordsys


// End of file 
