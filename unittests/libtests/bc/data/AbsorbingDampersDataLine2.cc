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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "AbsorbingDampersDataLine2.hh"

const char* pylith::bc::AbsorbingDampersDataLine2::_meshFilename = 
  "data/line2.mesh";

const int pylith::bc::AbsorbingDampersDataLine2::_numBasis = 1;
const int pylith::bc::AbsorbingDampersDataLine2::_numQuadPts = 1;
const double pylith::bc::AbsorbingDampersDataLine2::_quadPts[] = {
  0.0,
};
const double pylith::bc::AbsorbingDampersDataLine2::_quadWts[] = {
  1.0,
};
const double pylith::bc::AbsorbingDampersDataLine2::_basis[] = {
  1.0,
};
const double pylith::bc::AbsorbingDampersDataLine2::_basisDerivRef[] = {
  1.0,
};

const char* pylith::bc::AbsorbingDampersDataLine2::_spatialDBFilename = "data/elasticstrain1d.spatialdb";
const int pylith::bc::AbsorbingDampersDataLine2::_id = 2;
const char* pylith::bc::AbsorbingDampersDataLine2::_label = "bc0";

const double pylith::bc::AbsorbingDampersDataLine2::_dt =   0.25;
const double pylith::bc::AbsorbingDampersDataLine2::_fieldTmdt[] = {
  1.0,
  1.1,
  1.2,
};
const double pylith::bc::AbsorbingDampersDataLine2::_fieldT[] = {
  1.1,
  1.3,
  1.5,
};
const double pylith::bc::AbsorbingDampersDataLine2::_fieldTIncr[] = {
  1.2,
  1.5,
  1.8,
};

const int pylith::bc::AbsorbingDampersDataLine2::_spaceDim = 1;
const int pylith::bc::AbsorbingDampersDataLine2::_cellDim = 0;
const int pylith::bc::AbsorbingDampersDataLine2::_numVertices = 2;
const int pylith::bc::AbsorbingDampersDataLine2::_numCells = 2;
const int pylith::bc::AbsorbingDampersDataLine2::_numCorners = 1;
const int pylith::bc::AbsorbingDampersDataLine2::_cells[] = {
  2,
  4,
};


const double pylith::bc::AbsorbingDampersDataLine2::_dampingConsts[] = {
  12.5e+6,
  17.5e+6,
};
const double pylith::bc::AbsorbingDampersDataLine2::_valsResidual[] = {
  -12.5e+6*(1.1+1.2-1.0)/0.5,
  0.0,
  -17.5e+6*(1.5+1.8-1.2)/0.5,
};
const double pylith::bc::AbsorbingDampersDataLine2::_valsJacobian[] = {
  12.5e+6/0.5, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 17.5e+6/0.5,
};

pylith::bc::AbsorbingDampersDataLine2::AbsorbingDampersDataLine2(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);

  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  quadPts = const_cast<double*>(_quadPts);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDerivRef = const_cast<double*>(_basisDerivRef);

  spatialDBFilename = const_cast<char*>(_spatialDBFilename);
  id = _id;
  label = const_cast<char*>(_label);

  dt = _dt;
  fieldTIncr = const_cast<double*>(_fieldTIncr);
  fieldT = const_cast<double*>(_fieldT);
  fieldTmdt = const_cast<double*>(_fieldTmdt);

  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numCorners = _numCorners;
  cells = const_cast<int*>(_cells);

  dampingConsts = const_cast<double*>(_dampingConsts);
  valsResidual = const_cast<double*>(_valsResidual);
  valsJacobian = const_cast<double*>(_valsJacobian);
} // constructor

pylith::bc::AbsorbingDampersDataLine2::~AbsorbingDampersDataLine2(void)
{}


// End of file
