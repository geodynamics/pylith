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

#include "GeomDataHex3D.hh"

const int pylith::feassemble::GeomDataHex3D::_cellDim = 3;

const double pylith::feassemble::GeomDataHex3D::_gravityVec[] = {
  0.0, 0.0, -9.80665 };

const int pylith::feassemble::GeomDataHex3D::_spaceDim = 3;

const int pylith::feassemble::GeomDataHex3D::_numCorners = 8;

const int pylith::feassemble::GeomDataHex3D::_numLocs = 9;

const double pylith::feassemble::GeomDataHex3D::_vertices[] = {
  -1.3, -1.4, -0.8,
  1.2, -1.5, -0.9,
  1.4, 0.7, -1.2,
  -1.6, 0.4, -0.5,
  -1.7, -0.8, 1.8,
  2.1, -1.7, 0.6,
  2.3, 0.2, 1.9, 
  -1.8, 0.3, 2.2,
};

const double pylith::feassemble::GeomDataHex3D::_locations[] = {
 -1.0,-1.0,-1.0,
  1.0,-1.0,-1.0,
 -1.0, 1.0,-1.0,
  1.0, 1.0,-1.0,
 -1.0,-1.0, 1.0,
  1.0,-1.0, 1.0,
 -1.0, 1.0, 1.0,
  1.0, 1.0, 1.0,
 -0.6, 0.6, 0.4
};

// Reference cell has volume of 2.0, so divide by 2.0;
const double pylith::feassemble::GeomDataHex3D::_jacobian[] = {
  2.5/2.0, -0.3/2.0, -0.4/2.0, -0.1/2.0, 1.8/2.0, 0.6/2.0, -0.1/2.0, 0.3/2.0, 2.6/2.0,
  2.5/2.0, 0.2/2.0, 0.9/2.0, -0.1/2.0, 2.2/2.0, -0.2/2.0, -0.1/2.0, -0.3/2.0, 1.5/2.0,
  3.0/2.0, -0.3/2.0, -0.2/2.0, 0.3/2.0, 1.8/2.0, -0.1/2.0, -0.7/2.0, 0.3/2.0, 2.7/2.0,
  3.0/2.0, 0.2/2.0, 0.9/2.0, 0.3/2.0, 2.2/2.0, -0.5/2.0, -0.7/2.0, -0.3/2.0, 3.1/2.0,
  3.8/2.0, -0.1/2.0, -0.4/2.0, -0.9/2.0, 1.1/2.0, 0.6/2.0, -1.2/2.0, 0.4/2.0, 2.6/2.0,
  3.8/2.0, 0.2/2.0, 0.9/2.0, -0.9/2.0, 1.9/2.0, -0.2/2.0, -1.2/2.0, 1.3/2.0, 1.5/2.0,
  4.1/2.0, -0.1/2.0, -0.2/2.0, -0.1/2.0, 1.1/2.0, -0.1/2.0, -0.3/2.0, 0.4/2.0, 2.7/2.0,
  4.1/2.0, 0.2/2.0, 0.9/2.0, -0.1/2.0, 1.9/2.0, -0.5/2.0, -0.3/2.0, 1.3/2.0, 3.1/2.0,
  3.698/2.0, -0.088/2.0, -0.012/2.0, -0.116/2.0, 1.446/2.0, -0.056/2.0, -0.51/2.0, 0.46/2.0, 2.7/2.0
};

// Reference cell has volume of 2.0, so divide by 8.0;
const double pylith::feassemble::GeomDataHex3D::_jacobianDet[] = {
  11.13/8.0,
  8.359/8.0,
  14.622/8.0,
  21.199/8.0,
  9.41/8.0,
  13.135/8.0,
  12.253/8.0,
  27.302/8.0,
  14.494708/8.0
};

pylith::feassemble::GeomDataHex3D::GeomDataHex3D(void)
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

pylith::feassemble::GeomDataHex3D::~GeomDataHex3D(void)
{}


// End of file
