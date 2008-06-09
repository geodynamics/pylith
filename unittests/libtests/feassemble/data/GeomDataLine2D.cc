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

#include "GeomDataLine2D.hh"

const int pylith::feassemble::GeomDataLine2D::_cellDim = 1;

const int pylith::feassemble::GeomDataLine2D::_spaceDim = 2;

const double pylith::feassemble::GeomDataLine2D::_gravityVec[] = {
  0.0, -9.80665 };

const int pylith::feassemble::GeomDataLine2D::_numCorners = 2;

const int pylith::feassemble::GeomDataLine2D::_numLocs = 2;

const double pylith::feassemble::GeomDataLine2D::_vertices[] = {
  1.2, 2.4,
  4.5, -1.4
};

const double pylith::feassemble::GeomDataLine2D::_locations[] = {
  0.345,
  0.459
};

// Reference cell has length of 2.0, so divide by 2.0;
const double pylith::feassemble::GeomDataLine2D::_jacobian[] = {
  3.3/2.0, -3.8/2.0,
  3.3/2.0, -3.8/2.0
};

// Reference cell has length of 2.0, so divide by 2.0;
const double pylith::feassemble::GeomDataLine2D::_jacobianDet[] = {
  5.0328918128646478/2.0,
  5.0328918128646478/2.0
};

pylith::feassemble::GeomDataLine2D::GeomDataLine2D(void)
{ // constructor
  cellDim = _cellDim;
  spaceDim = _spaceDim;
  numCorners = _numCorners;
  numLocs = _numLocs;
  gravityVec = const_cast<double*>(_gravityVec);
  vertices = const_cast<double*>(_vertices);
  locations = const_cast<double*>(_locations);
  jacobian = const_cast<double*>(_jacobian);
  jacobianDet = const_cast<double*>(_jacobianDet);
} // constructor

pylith::feassemble::GeomDataLine2D::~GeomDataLine2D(void)
{}


// End of file
