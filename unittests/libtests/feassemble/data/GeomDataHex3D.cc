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

const double pylith::feassemble::GeomDataHex3D::_jacobian[] = {
  2.5, -0.3, -0.4, -0.1, 1.8, 0.6, -0.1, 0.3, 2.6,
  2.5, 0.2, 0.9, -0.1, 2.2, -0.2, -0.1, -0.3, 1.5,
  3.0, -0.3, -0.2, 0.3, 1.8, -0.1, -0.7, 0.3, 2.7,
  3.0, 0.2, 0.9, 0.3, 2.2, -0.5, -0.7, -0.3, 3.1,
  3.8, -0.1, -0.4, -0.9, 1.1, 0.6, -1.2, 0.4, 2.6,
  3.8, 0.2, 0.9, -0.9, 1.9, -0.2, -1.2, 1.3, 1.5,
  4.1, -0.1, -0.2, -0.1, 1.1, -0.1, -0.3, 0.4, 2.7,
  4.1, 0.2, 0.9, -0.1, 1.9, -0.5, -0.3, 1.3, 3.1,
  3.698, -0.088, -0.012, -0.116, 1.446, -0.056, -0.51, 0.46, 2.7
};

const double pylith::feassemble::GeomDataHex3D::_jacobianDet[] = {
  11.13,
  8.359,
  14.622,
  21.199,
  9.41,
  13.135,
  12.253,
  27.302,
  14.494708
};

pylith::feassemble::GeomDataHex3D::GeomDataHex3D(void)
{ // constructor
  cellDim = _cellDim;
  spaceDim = _spaceDim;
  numCorners = _numCorners;
  numLocs = _numLocs;
  vertices = const_cast<double*>(_vertices);
  locations = const_cast<double*>(_locations);
  jacobian = const_cast<double*>(_jacobian);
  jacobianDet = const_cast<double*>(_jacobianDet);
} // constructor

pylith::feassemble::GeomDataHex3D::~GeomDataHex3D(void)
{}


// End of file
