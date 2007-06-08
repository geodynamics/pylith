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

#include "GeomDataTri2D.hh"

const int pylith::feassemble::GeomDataTri2D::_cellDim = 2;

const int pylith::feassemble::GeomDataTri2D::_spaceDim = 2;

const int pylith::feassemble::GeomDataTri2D::_numCorners = 3;

const int pylith::feassemble::GeomDataTri2D::_numLocs = 2;

const double pylith::feassemble::GeomDataTri2D::_vertices[] = {
  1.0,  1.2,
  3.0,  2.0,
  1.5,  4.0
};

const double pylith::feassemble::GeomDataTri2D::_locations[] = {
  0.345, 0.397,
  0.459, 0.727
};

const double pylith::feassemble::GeomDataTri2D::_jacobian[] = {
  2.0, 0.5, 0.8, 2.8,
  2.0, 0.5, 0.8, 2.8
};

const double pylith::feassemble::GeomDataTri2D::_jacobianDet[] = {
  5.2,
  5.2
};

pylith::feassemble::GeomDataTri2D::GeomDataTri2D(void)
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

pylith::feassemble::GeomDataTri2D::~GeomDataTri2D(void)
{}


// End of file
