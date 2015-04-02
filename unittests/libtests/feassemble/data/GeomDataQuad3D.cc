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

#include "GeomDataQuad3D.hh"

const int pylith::feassemble::GeomDataQuad3D::_cellDim = 2;

const int pylith::feassemble::GeomDataQuad3D::_spaceDim = 3;

const PylithScalar pylith::feassemble::GeomDataQuad3D::_gravityVec[] = {
  0.0, 0.0, -9.80665 };

const int pylith::feassemble::GeomDataQuad3D::_numCorners = 4;

const int pylith::feassemble::GeomDataQuad3D::_numLocs = 5;

const PylithScalar pylith::feassemble::GeomDataQuad3D::_vertices[] = {
  1.1, -0.9, 0.1,
  1.0, 0.7, 0.8,
  -0.4, 0.6, 1.2,
  -0.3, -1.3, -0.1,
};

const PylithScalar pylith::feassemble::GeomDataQuad3D::_locations[] = {
 -1.0, -1.0,
  1.0, -1.0,
 -1.0,  1.0,
  1.0,  1.0,
 -0.4,  0.6
};

// Reference cell has area of 2.0, so divide by 2.0;
const PylithScalar pylith::feassemble::GeomDataQuad3D::_jacobian[] = {
  -0.1/2.0, -1.4/2.0, 1.6/2.0, -0.4/2.0, 0.7/2.0, -0.2/2.0,
  -0.1/2.0, -1.4/2.0, 1.6/2.0, -0.1/2.0, 0.7/2.0, 0.4/2.0,
  -0.1/2.0, -1.4/2.0, 1.9/2.0, -0.4/2.0, 1.3/2.0, -0.2/2.0,
  -0.1/2.0, -1.4/2.0, 1.9/2.0, -0.1/2.0, 1.3/2.0, 0.4/2.0,
  -0.1/2.0, -1.4/2.0, 1.84/2.0, -0.31/2.0, 1.18/2.0, -0.02/2.0
};

// Reference cell has area of 2.0, so divide by 4.0;
const PylithScalar pylith::feassemble::GeomDataQuad3D::_jacobianDet[] = {
  2.4899799195977463/4.0,
  2.5397243944963788/4.0,
  3.2703516630478742/4.0,
  3.3300750742288079/4.0,
  3.1049003204611898/4.0
};

pylith::feassemble::GeomDataQuad3D::GeomDataQuad3D(void)
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

pylith::feassemble::GeomDataQuad3D::~GeomDataQuad3D(void)
{}


// End of file
