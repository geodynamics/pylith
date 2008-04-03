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

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriter::DataWriter(void) :
  _numTimeSteps(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriter::~DataWriter(void)
{ // destructor
} // destructor  

// ----------------------------------------------------------------------
// Prepare for writing files.
void
pylith::meshio::DataWriter::open(
			       const ALE::Obj<Mesh>& mesh,
			       const spatialdata::geocoords::CoordSys* csMesh,
			       const int numTimeSteps,
			       const char* label,
			       const int labelId)
{ // open
  _numTimeSteps = numTimeSteps;
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
			       const ALE::Obj<Mesh>& mesh,
			       const spatialdata::geocoords::CoordSys* csMesh,
			       const char* label,
			       const int labelId)
{ // openTimeStep
} // openTimeStep

// ----------------------------------------------------------------------
// Cleanup after writing data for a time step.
void
pylith::meshio::DataWriter::closeTimeStep(void)
{ // closeTimeStep
} // closeTimeStep

// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::DataWriter::DataWriter(const DataWriter& w)
{ // copy constructor
} // copy constructor


// End of file 
