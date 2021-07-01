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

#include "AbsorbingDampersDataTet4.hh"

const char* pylith::bc::AbsorbingDampersDataTet4::_meshFilename = 
  "data/tet4.mesh";

const int pylith::bc::AbsorbingDampersDataTet4::_numBasis = 3;
const int pylith::bc::AbsorbingDampersDataTet4::_numQuadPts = 1;
const PylithScalar pylith::bc::AbsorbingDampersDataTet4::_quadPts[] = {
  -0.3333333333333333, -0.3333333333333333
};
const PylithScalar pylith::bc::AbsorbingDampersDataTet4::_quadWts[] = {
  2.0,
};
const PylithScalar pylith::bc::AbsorbingDampersDataTet4::_basis[] = {
  0.3333333333333333,
  0.3333333333333333,
  0.3333333333333333,
};
const PylithScalar pylith::bc::AbsorbingDampersDataTet4::_basisDerivRef[] = {
 -0.5, -0.5,
  0.5,  0.0,
  0.0,  0.5,
};

const char* pylith::bc::AbsorbingDampersDataTet4::_spatialDBFilename = 
  "data/elasticisotropic3d.spatialdb";
const int pylith::bc::AbsorbingDampersDataTet4::_id = 2;
const char* pylith::bc::AbsorbingDampersDataTet4::_label = "bc2";

const PylithScalar pylith::bc::AbsorbingDampersDataTet4::_dt =   0.25;
const PylithScalar pylith::bc::AbsorbingDampersDataTet4::_fieldTmdt[] = {
  1.0,  2.4,  3.0,
  1.1,  1.8,  3.2,
  1.2,  2.4,  3.4,
  1.3,  2.2,  3.6
};
const PylithScalar pylith::bc::AbsorbingDampersDataTet4::_fieldT[] = {
  1.1,  2.0,  3.2,
  1.3,  2.1,  3.6,
  1.5,  2.2,  4.0,
  1.7,  2.3,  4.4,
};
const PylithScalar pylith::bc::AbsorbingDampersDataTet4::_fieldTIncr[] = {
  1.2,  1.6,  3.4,
  1.5,  2.4,  4.0,
  1.8,  2.0,  4.6,
  2.1,  2.4,  5.2
};

const int pylith::bc::AbsorbingDampersDataTet4::_spaceDim = 3;
const int pylith::bc::AbsorbingDampersDataTet4::_cellDim = 2;
const int pylith::bc::AbsorbingDampersDataTet4::_numVertices = 3;
const int pylith::bc::AbsorbingDampersDataTet4::_numCells = 1;
const int pylith::bc::AbsorbingDampersDataTet4::_numCorners = 3;
/* Now vertices are renumbered in the submesh */
const int pylith::bc::AbsorbingDampersDataTet4::_cells[] = {
  3 /*5*/,  2 /*4*/,  1 /*2*/,
};


const PylithScalar pylith::bc::AbsorbingDampersDataTet4::_dampingConsts[] = {
  1.25e+07,  7.5e+06,  7.5e+06
};
const PylithScalar pylith::bc::AbsorbingDampersDataTet4::_valsResidual[] = {
  -8.19444444e+06, -4.58333333e+06, -1.23333333e+07,
  0.0,              0.0,             0.0,
  -8.19444444e+06, -4.58333333e+06, -1.23333333e+07,
  -8.19444444e+06, -4.58333333e+06, -1.23333333e+07,
  0.0,              0.0,             0.0,  
};
const PylithScalar pylith::bc::AbsorbingDampersDataTet4::_valsJacobian[] = {
  1.38888889e+06, 0.0, 0.0, // 0x
  0.0, 0.0, 0.0,
  1.38888889e+06, 0.0, 0.0,
  1.38888889e+06, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 8.33333333e+05, 0.0, // 0y
  0.0, 0.0, 0.0,
  0.0, 8.33333333e+05, 0.0,
  0.0, 8.33333333e+05, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 8.33333333e+05, // 0z
  0.0, 0.0, 0.0,
  0.0, 0.0, 8.33333333e+05,
  0.0, 0.0, 8.33333333e+05,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 1x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 1y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 1z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  1.38888889e+06, 0.0, 0.0, // 2x
  0.0, 0.0, 0.0,
  1.38888889e+06, 0.0, 0.0,
  1.38888889e+06, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 8.33333333e+05, 0.0, // 2y
  0.0, 0.0, 0.0,
  0.0, 8.33333333e+05, 0.0,
  0.0, 8.33333333e+05, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 8.33333333e+05, // 2z
  0.0, 0.0, 0.0,
  0.0, 0.0, 8.33333333e+05,
  0.0, 0.0, 8.33333333e+05,
  0.0, 0.0, 0.0,
  1.38888889e+06, 0.0, 0.0, // 3x
  0.0, 0.0, 0.0,
  1.38888889e+06, 0.0, 0.0,
  1.38888889e+06, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 8.33333333e+05, 0.0, // 3y
  0.0, 0.0, 0.0,
  0.0, 8.33333333e+05, 0.0,
  0.0, 8.33333333e+05, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 8.33333333e+05, // 3z
  0.0, 0.0, 0.0,
  0.0, 0.0, 8.33333333e+05,
  0.0, 0.0, 8.33333333e+05,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 4x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 4y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 4z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
};

pylith::bc::AbsorbingDampersDataTet4::AbsorbingDampersDataTet4(void)
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

pylith::bc::AbsorbingDampersDataTet4::~AbsorbingDampersDataTet4(void)
{}


// End of file
