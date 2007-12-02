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

#include "GeomDataLine3D.hh"

const int pylith::feassemble::GeomDataLine3D::_cellDim = 1;

const int pylith::feassemble::GeomDataLine3D::_spaceDim = 3;

const int pylith::feassemble::GeomDataLine3D::_numCorners = 2;

const int pylith::feassemble::GeomDataLine3D::_numLocs = 2;

const double pylith::feassemble::GeomDataLine3D::_vertices[] = {
  1.2, 2.3, 2.4,
  4.5, -1.5, 5.3
};

const double pylith::feassemble::GeomDataLine3D::_locations[] = {
  0.345,
  0.459
};

// Reference cell has length of 2.0, so divide by 2.0;
const double pylith::feassemble::GeomDataLine3D::_jacobian[] = {
  3.3/2.0, -3.8/2.0, 2.9/2.0,
  3.3/2.0, -3.8/2.0, 2.9/2.0
};

// Reference cell has length of 2.0, so divide by 2.0;
const double pylith::feassemble::GeomDataLine3D::_jacobianDet[] = {
  5.8086142925830417/2.0,
  5.8086142925830417/2.0
};

pylith::feassemble::GeomDataLine3D::GeomDataLine3D(void)
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

pylith::feassemble::GeomDataLine3D::~GeomDataLine3D(void)
{}


// End of file
