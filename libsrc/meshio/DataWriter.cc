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
pylith::meshio::DataWriter::DataWriter(void)
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
