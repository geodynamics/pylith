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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "GeomDataLine1D.hh"

const int pylith::feassemble::GeomDataLine1D::_cellDim = 1;

const int pylith::feassemble::GeomDataLine1D::_spaceDim = 1;

const PylithScalar pylith::feassemble::GeomDataLine1D::_gravityVec[] = {
  -9.80665 };

const int pylith::feassemble::GeomDataLine1D::_numCorners = 2;

const int pylith::feassemble::GeomDataLine1D::_numLocs = 2;

const PylithScalar pylith::feassemble::GeomDataLine1D::_vertices[] = {
  1.2,
  4.5
};

const PylithScalar pylith::feassemble::GeomDataLine1D::_locations[] = {
  0.345,
  0.459
};

// Reference cell has length of 2.0, so divide by 2.0;
const PylithScalar pylith::feassemble::GeomDataLine1D::_jacobian[] = {
  3.3/2.0,
  3.3/2.0
};

// Reference cell has length of 2.0, so divide by 2.0;
const PylithScalar pylith::feassemble::GeomDataLine1D::_jacobianDet[] = {
  3.3/2.0,
  3.3/2.0
};

pylith::feassemble::GeomDataLine1D::GeomDataLine1D(void)
{ // constructor
  cellDim = _cellDim;
  spaceDim = _spaceDim;
  numCorners = _numCorners;
  numLocs = _numLocs;
  gravityVec = const_cast<PylithScalar*>(_gravityVec);
  vertices = const_cast<PylithScalar*>(_vertices);
  locations = const_cast<PylithScalar*>(_locations);
  jacobian = const_cast<PylithScalar*>(_jacobian);
  jacobianDet = const_cast<PylithScalar*>(_jacobianDet);
} // constructor

pylith::feassemble::GeomDataLine1D::~GeomDataLine1D(void)
{}


// End of file
