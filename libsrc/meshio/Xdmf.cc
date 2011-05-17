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

#include "Xdmf.hh" // implementation of class methods

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::meshio::Xdmf::Xdmf(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::Xdmf::~Xdmf(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Write Xdmf file associated with HDF5 file.
void
pylith::meshio::Xdmf::write(const char* filenameXdfm,
			    const char* filenameHDF5)
{ // write
} // write

// ----------------------------------------------------------------------
// Write domain cell information.
void
pylith::meshio::Xdmf::_writeDomainCells(const int numcells,
					const int numCorners)
{ // _writeDomainCells
} // _writeDomainCells

// ----------------------------------------------------------------------
// Write domain vertices information.
void
pylith::meshio::Xdmf::_writeDomainVertices(const int numVertices,
					   const int spaceDim)
{ // _writeDomainVertices
} // _writeDomainVertices

// ----------------------------------------------------------------------
// Write grid topology information.
void
pylith::meshio::Xdmf::_writeGridTopology(const char* cellType,
					 const int numells)
{ // _writeGridTopology
} // _writeGridTopology

// ----------------------------------------------------------------------
// Write Grid geometry.
void
pylith::meshio::Xdmf::_writeGridGeometry(const int spaceDim)
{ // _writeGridGeometry
} // _writeGridGeometry

// ----------------------------------------------------------------------
// Write grid attribute.
void
pylith::meshio::Xdmf::_writeGridAttribute(const char* name,
					  const char* center,
					  const int numTimeSteps,
					  const int numPoints,
					  const int fiberDim,
					  const int iTime)
{ // _writeGridAttribute
} // _writeGridAttribute


// End of file
