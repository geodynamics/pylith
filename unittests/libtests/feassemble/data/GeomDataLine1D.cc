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

#include "GeomDataLine1D.hh"

const int pylith::feassemble::GeomDataLine1D::_cellDim = 1;

const int pylith::feassemble::GeomDataLine1D::_spaceDim = 1;

const int pylith::feassemble::GeomDataLine1D::_numCorners = 2;

const int pylith::feassemble::GeomDataLine1D::_numLocs = 2;

const double pylith::feassemble::GeomDataLine1D::_vertices[] = {
  1.2,
  4.5
};

const double pylith::feassemble::GeomDataLine1D::_locations[] = {
  0.345,
  0.459
};

const double pylith::feassemble::GeomDataLine1D::_jacobian[] = {
  3.3,
  3.3
};

const double pylith::feassemble::GeomDataLine1D::_jacobianDet[] = {
  3.3,
  3.3
};

pylith::feassemble::GeomDataLine1D::GeomDataLine1D(void)
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

pylith::feassemble::GeomDataLine1D::~GeomDataLine1D(void)
{}


// End of file
