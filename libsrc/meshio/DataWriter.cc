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

#include "DataWriter.hh" // implementation of class methods

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriter::DataWriter(void) :
  _cs(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriter::~DataWriter(void)
{ // destructor
  delete _cs; _cs = 0;
} // destructor  

// ----------------------------------------------------------------------
// Set coordinate system for output.
void
pylith::meshio::DataWriter::coordsys(const spatialdata::geocoords::CoordSys* cs)
{ // coordsys
  delete _cs; _cs = (0 != cs) ? cs->clone() : 0;
} // coordsys

// ----------------------------------------------------------------------
// Prepare for writing files.
void
pylith::meshio::DataWriter::open(
			       const ALE::Obj<ALE::Mesh>& mesh,
			       const spatialdata::geocoords::CoordSys* csMesh)
{ // open
} // open

// ----------------------------------------------------------------------
// Close output files.
void
pylith::meshio::DataWriter::close(void)
{ // close
} // close

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
void
pylith::meshio::DataWriter::openTimeStep(
			       const double t,
			       const ALE::Obj<ALE::Mesh>& mesh,
			       const spatialdata::geocoords::CoordSys* csMesh)
{ // openTimeStep
} // openTimeStep

// ----------------------------------------------------------------------
/// Cleanup after writing data for a time step.
void
pylith::meshio::DataWriter::closeTimeStep(void)
{ // closeTimeStep
} // closeTimeStep


// End of file 
