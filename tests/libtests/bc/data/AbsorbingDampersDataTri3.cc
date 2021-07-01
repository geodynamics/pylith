// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include "AbsorbingDampersDataTri3.hh"

const char* pylith::bc::AbsorbingDampersDataTri3::_meshFilename = 
  "data/tri3.mesh";

const int pylith::bc::AbsorbingDampersDataTri3::_numBasis = 2;
const int pylith::bc::AbsorbingDampersDataTri3::_numQuadPts = 1;
const PylithScalar pylith::bc::AbsorbingDampersDataTri3::_quadPts[] = {
  0.0,
};
const PylithScalar pylith::bc::AbsorbingDampersDataTri3::_quadWts[] = {
  2.0,
};
const PylithScalar pylith::bc::AbsorbingDampersDataTri3::_basis[] = {
  0.5,
  0.5,
};
const PylithScalar pylith::bc::AbsorbingDampersDataTri3::_basisDerivRef[] = {
  -0.5,
   0.5,
};

const char* pylith::bc::AbsorbingDampersDataTri3::_spatialDBFilename = 
  "data/elasticplanestrain.spatialdb";
const int pylith::bc::AbsorbingDampersDataTri3::_id = 2;
const char* pylith::bc::AbsorbingDampersDataTri3::_label = "bc";

const PylithScalar pylith::bc::AbsorbingDampersDataTri3::_dt =   0.25;
const PylithScalar pylith::bc::AbsorbingDampersDataTri3::_fieldTmdt[] = {
  1.0,  2.4,
  1.1,  1.8,
  1.2,  2.4,
  1.3,  2.2,
};
const PylithScalar pylith::bc::AbsorbingDampersDataTri3::_fieldT[] = {
  1.1,  2.0,
  1.3,  2.1,
  1.5,  2.2,
  1.7,  2.3,
};
const PylithScalar pylith::bc::AbsorbingDampersDataTri3::_fieldTIncr[] = {
  1.2,  1.6,
  1.5,  2.4,
  1.8,  2.0,
  2.1,  2.4,
};

const int pylith::bc::AbsorbingDampersDataTri3::_spaceDim = 2;
const int pylith::bc::AbsorbingDampersDataTri3::_cellDim = 1;
const int pylith::bc::AbsorbingDampersDataTri3::_numVertices = 2;
const int pylith::bc::AbsorbingDampersDataTri3::_numCells = 1;
const int pylith::bc::AbsorbingDampersDataTri3::_numCorners = 2;
/* Now vertices are renumbered in the submesh */
const int pylith::bc::AbsorbingDampersDataTri3::_cells[] = {
  1 /*3*/,
  2 /*5*/,
};


const PylithScalar pylith::bc::AbsorbingDampersDataTri3::_dampingConsts[] = {
  1.41421356e+07,  3.53553391e+06
};
const PylithScalar pylith::bc::AbsorbingDampersDataTri3::_valsResidual[] = {
  0.0, 0.0,
  -4.20000000e+07,   -1.30000000e+07,
  0.0, 0.0,
  -4.20000000e+07,   -1.30000000e+07,
};
const PylithScalar pylith::bc::AbsorbingDampersDataTri3::_valsJacobian[] = {
  0.0, 0.0, // 0x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 0y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 1x
  1.0e+07, 0.0,
  0.0, 0.0,
  1.0e+07, 0.0,
  0.0, 0.0, // 1y
  0.0, 2.5e+06,
  0.0, 0.0,
  0.0, 2.5e+06,
  0.0, 0.0, // 2x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 2y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 3x
  1.0e+07, 0.0,
  0.0, 0.0,
  1.0e+07, 0.0,
  0.0, 0.0, // 3y
  0.0, 2.5e+06,
  0.0, 0.0,
  0.0, 2.5e+06,
};

pylith::bc::AbsorbingDampersDataTri3::AbsorbingDampersDataTri3(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);

  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDerivRef = const_cast<PylithScalar*>(_basisDerivRef);

  spatialDBFilename = const_cast<char*>(_spatialDBFilename);
  id = _id;
  label = const_cast<char*>(_label);

  dt = _dt;
  fieldTIncr = const_cast<PylithScalar*>(_fieldTIncr);
  fieldT = const_cast<PylithScalar*>(_fieldT);
  fieldTmdt = const_cast<PylithScalar*>(_fieldTmdt);

  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numCorners = _numCorners;
  cells = const_cast<int*>(_cells);

  dampingConsts = const_cast<PylithScalar*>(_dampingConsts);
  valsResidual = const_cast<PylithScalar*>(_valsResidual);
  valsJacobian = const_cast<PylithScalar*>(_valsJacobian);
} // constructor

pylith::bc::AbsorbingDampersDataTri3::~AbsorbingDampersDataTri3(void)
{}


// End of file
