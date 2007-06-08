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

#include "GeomDataPoint2D.hh"

const int pylith::feassemble::GeomDataPoint2D::_cellDim = 0;

const int pylith::feassemble::GeomDataPoint2D::_spaceDim = 2;

const int pylith::feassemble::GeomDataPoint2D::_numCorners = 1;

const int pylith::feassemble::GeomDataPoint2D::_numLocs = 2;

const double pylith::feassemble::GeomDataPoint2D::_vertices[] = {
  1.3, 5.4,
  4.1, 7.5
};

const double pylith::feassemble::GeomDataPoint2D::_locations[] = {
  0.0,
  0.0
};

const double pylith::feassemble::GeomDataPoint2D::_jacobian[] = {
  1.0,
  1.0
};

const double pylith::feassemble::GeomDataPoint2D::_jacobianDet[] = {
  1.0,
  1.0
};

pylith::feassemble::GeomDataPoint2D::GeomDataPoint2D(void)
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

pylith::feassemble::GeomDataPoint2D::~GeomDataPoint2D(void)
{}


// End of file
