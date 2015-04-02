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

#include "GeomDataQuad2D.hh"

const int pylith::feassemble::GeomDataQuad2D::_cellDim = 2;

const int pylith::feassemble::GeomDataQuad2D::_spaceDim = 2;

const PylithScalar pylith::feassemble::GeomDataQuad2D::_gravityVec[] = {
  0.0, -9.80665 };

const int pylith::feassemble::GeomDataQuad2D::_numCorners = 4;

const int pylith::feassemble::GeomDataQuad2D::_numLocs = 5;

const PylithScalar pylith::feassemble::GeomDataQuad2D::_vertices[] = {
  0.3, 0.1,
  0.8, -0.2,
  0.7, 1.2,
  -0.1, 1.6,
};

const PylithScalar pylith::feassemble::GeomDataQuad2D::_locations[] = {
 -1.0, -1.0,
  1.0, -1.0,
 -1.0,  1.0,
  1.0,  1.0,
 -0.2,  0.4
};

// Reference cell has area of 2.0, so divide by 2.0;
const PylithScalar pylith::feassemble::GeomDataQuad2D::_jacobian[] = {
  0.5/2.0, -0.4/2.0, -0.3/2.0, 1.5/2.0,
  0.5/2.0, -0.1/2.0, -0.3/2.0, 1.4/2.0,
  0.8/2.0, -0.4/2.0, -0.4/2.0, 1.5/2.0,
  0.8/2.0, -0.1/2.0, -0.4/2.0, 1.4/2.0,
  0.71/2.0, -0.28/2.0, -0.37/2.0, 1.46/2.0
};

// Reference cell has area of 2.0, so divide by 4.0;
const PylithScalar pylith::feassemble::GeomDataQuad2D::_jacobianDet[] = {
  0.63/4.0,
  0.67/4.0,
  1.04/4.0,
  1.08/4.0,
  0.933/4.0
};

pylith::feassemble::GeomDataQuad2D::GeomDataQuad2D(void)
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

pylith::feassemble::GeomDataQuad2D::~GeomDataQuad2D(void)
{}


// End of file
