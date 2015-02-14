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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "GeomDataLine3D.hh"

const int pylith::feassemble::GeomDataLine3D::_cellDim = 1;

const int pylith::feassemble::GeomDataLine3D::_spaceDim = 3;

const PylithScalar pylith::feassemble::GeomDataLine3D::_gravityVec[] = {
  0.0, 0.0, -9.80665 };

const int pylith::feassemble::GeomDataLine3D::_numCorners = 2;

const int pylith::feassemble::GeomDataLine3D::_numLocs = 2;

const PylithScalar pylith::feassemble::GeomDataLine3D::_vertices[] = {
  1.2, 2.3, 2.4,
  4.5, -1.5, 5.3
};

const PylithScalar pylith::feassemble::GeomDataLine3D::_locations[] = {
  0.345,
  0.459
};

// Reference cell has length of 2.0, so divide by 2.0;
const PylithScalar pylith::feassemble::GeomDataLine3D::_jacobian[] = {
  3.3/2.0, -3.8/2.0, 2.9/2.0,
  3.3/2.0, -3.8/2.0, 2.9/2.0
};

// Reference cell has length of 2.0, so divide by 2.0;
const PylithScalar pylith::feassemble::GeomDataLine3D::_jacobianDet[] = {
  5.8086142925830417/2.0,
  5.8086142925830417/2.0
};

pylith::feassemble::GeomDataLine3D::GeomDataLine3D(void)
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

pylith::feassemble::GeomDataLine3D::~GeomDataLine3D(void)
{}


// End of file
