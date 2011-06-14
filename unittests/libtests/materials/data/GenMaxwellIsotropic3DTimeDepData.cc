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

// DO NOT EDIT THIS FILE
// This file was generated from python application genmaxwellisotropic3dtimedep.

#include "GenMaxwellIsotropic3DTimeDepData.hh"

const int pylith::materials::GenMaxwellIsotropic3DTimeDepData::_dimension = 3;

const int pylith::materials::GenMaxwellIsotropic3DTimeDepData::_numLocs = 2;

const int pylith::materials::GenMaxwellIsotropic3DTimeDepData::_numProperties = 9;

const int pylith::materials::GenMaxwellIsotropic3DTimeDepData::_numStateVars = 4;

const int pylith::materials::GenMaxwellIsotropic3DTimeDepData::_numDBProperties = 9;

const int pylith::materials::GenMaxwellIsotropic3DTimeDepData::_numDBStateVars = 24;

const int pylith::materials::GenMaxwellIsotropic3DTimeDepData::_numPropsQuadPt = 9;

const int pylith::materials::GenMaxwellIsotropic3DTimeDepData::_numVarsQuadPt = 24;

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_lengthScale =   1.00000000e+03;

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_timeScale =   1.00000000e+00;

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_pressureScale =   2.25000000e+10;

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_densityScale =   1.00000000e+03;

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_dtStableImplicit =   8.88888889e+06;

const int pylith::materials::GenMaxwellIsotropic3DTimeDepData::_numPropertyValues[] = {
1,
1,
1,
1,
1,
1,
1,
1,
1,
};

const int pylith::materials::GenMaxwellIsotropic3DTimeDepData::_numStateVarValues[] = {
6,
6,
6,
6,
};

const char* pylith::materials::GenMaxwellIsotropic3DTimeDepData::_dbPropertyValues[] = {
"density",
"vs",
"vp",
"shear-ratio-1",
"shear-ratio-2",
"shear-ratio-3",
"viscosity-1",
"viscosity-2",
"viscosity-3",
};

const char* pylith::materials::GenMaxwellIsotropic3DTimeDepData::_dbStateVarValues[] = {
"total-strain-xx",
"total-strain-yy",
"total-strain-zz",
"total-strain-xy",
"total-strain-yz",
"total-strain-xz",
"viscous-strain-1-xx",
"viscous-strain-1-yy",
"viscous-strain-1-zz",
"viscous-strain-1-xy",
"viscous-strain-1-yz",
"viscous-strain-1-xz",
"viscous-strain-2-xx",
"viscous-strain-2-yy",
"viscous-strain-2-zz",
"viscous-strain-2-xy",
"viscous-strain-2-yz",
"viscous-strain-2-xz",
"viscous-strain-3-xx",
"viscous-strain-3-yy",
"viscous-strain-3-zz",
"viscous-strain-3-xy",
"viscous-strain-3-yz",
"viscous-strain-3-xz",
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_dbProperties[] = {
  2.50000000e+03,
  3.00000000e+03,
  5.19615242e+03,
  5.00000000e-01,
  1.00000000e-01,
  2.00000000e-01,
  1.00000000e+18,
  1.00000000e+17,
  1.00000000e+19,
  2.00000000e+03,
  1.20000000e+03,
  2.07846097e+03,
  2.00000000e-01,
  2.00000000e-01,
  2.00000000e-01,
  1.00000000e+18,
  1.00000000e+19,
  1.00000000e+20,
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_dbStateVars[] = {
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_properties[] = {
  2.50000000e+03,
  2.25000000e+10,
  2.25000000e+10,
  5.00000000e-01,
  1.00000000e-01,
  2.00000000e-01,
  8.88888889e+07,
  4.44444444e+07,
  2.22222222e+09,
  2.00000000e+03,
  2.88000000e+09,
  2.88000000e+09,
  2.00000000e-01,
  2.00000000e-01,
  2.00000000e-01,
  1.73611111e+09,
  1.73611111e+10,
  1.73611111e+11,
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_stateVars[] = {
  1.10000000e-04,
  2.20000000e-04,
  3.30000000e-04,
  4.40000000e-04,
  5.50000000e-04,
  6.60000000e-04,
 -1.10000000e-04,
  0.00000000e+00,
  1.10000000e-04,
  4.40000000e-04,
  5.50000000e-04,
  6.60000000e-04,
 -1.10000000e-04,
  0.00000000e+00,
  1.10000000e-04,
  4.40000000e-04,
  5.50000000e-04,
  6.60000000e-04,
 -1.10000000e-04,
  0.00000000e+00,
  1.10000000e-04,
  4.40000000e-04,
  5.50000000e-04,
  6.60000000e-04,
  1.20000000e-04,
  2.30000000e-04,
  3.40000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
 -1.10000000e-04,
 -2.71050543e-20,
  1.10000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
 -1.10000000e-04,
 -2.71050543e-20,
  1.10000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
 -1.10000000e-04,
 -2.71050543e-20,
  1.10000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_propertiesNondim[] = {
  2.50000000e+00,
  1.00000000e+00,
  1.00000000e+00,
  5.00000000e-01,
  1.00000000e-01,
  2.00000000e-01,
  8.88888889e+07,
  4.44444444e+07,
  2.22222222e+09,
  2.00000000e+00,
  1.28000000e-01,
  1.28000000e-01,
  2.00000000e-01,
  2.00000000e-01,
  2.00000000e-01,
  1.73611111e+09,
  1.73611111e+10,
  1.73611111e+11,
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_stateVarsNondim[] = {
  1.10000000e-04,
  2.20000000e-04,
  3.30000000e-04,
  4.40000000e-04,
  5.50000000e-04,
  6.60000000e-04,
 -1.10000000e-04,
  0.00000000e+00,
  1.10000000e-04,
  4.40000000e-04,
  5.50000000e-04,
  6.60000000e-04,
 -1.10000000e-04,
  0.00000000e+00,
  1.10000000e-04,
  4.40000000e-04,
  5.50000000e-04,
  6.60000000e-04,
 -1.10000000e-04,
  0.00000000e+00,
  1.10000000e-04,
  4.40000000e-04,
  5.50000000e-04,
  6.60000000e-04,
  1.20000000e-04,
  2.30000000e-04,
  3.40000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
 -1.10000000e-04,
 -2.71050543e-20,
  1.10000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
 -1.10000000e-04,
 -2.71050543e-20,
  1.10000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
 -1.10000000e-04,
 -2.71050543e-20,
  1.10000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_density[] = {
  2.50000000e+03,
  2.00000000e+03,
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_strain[] = {
  1.20000000e-04,
  2.30000000e-04,
  3.40000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
  1.30000000e-04,
  2.40000000e-04,
  3.50000000e-04,
  4.60000000e-04,
  5.70000000e-04,
  6.80000000e-04,
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_stress[] = {
  1.73848741e+07,
  2.23190000e+07,
  2.72531259e+07,
  1.99361456e+07,
  2.48702715e+07,
  2.98043974e+07,
  2.03492020e+06,
  2.66720000e+06,
  3.29947980e+06,
  2.55607698e+06,
  3.18835678e+06,
  3.82063657e+06,
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_elasticConsts[] = {
  6.74761292e+10,
  2.25119367e+10,
  2.25119367e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.25119367e+10,
  6.74761292e+10,
  2.25119367e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.25119386e+10,
  2.25119386e+10,
  6.74761292e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.49641924e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.49641924e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.49641906e+10,
  8.63995089e+09,
  2.88002449e+09,
  2.88002449e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.88002449e+09,
  8.63995100e+09,
  2.88002449e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.88002472e+09,
  2.88002472e+09,
  8.63995077e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.75992605e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.75992628e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.75992605e+09,
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_initialStress[] = {
  2.10000000e+04,
  2.20000000e+04,
  2.30000000e+04,
  2.40000000e+04,
  2.50000000e+04,
  2.60000000e+04,
  5.10000000e+04,
  5.20000000e+04,
  5.30000000e+04,
  5.40000000e+04,
  5.50000000e+04,
  5.60000000e+04,
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_initialStrain[] = {
  3.10000000e-05,
  3.20000000e-05,
  3.30000000e-05,
  3.40000000e-05,
  3.50000000e-05,
  3.60000000e-05,
  6.10000000e-05,
  6.20000000e-05,
  6.30000000e-05,
  6.40000000e-05,
  6.50000000e-05,
  6.60000000e-05,
};

const double pylith::materials::GenMaxwellIsotropic3DTimeDepData::_stateVarsUpdated[] = {
  1.20000000e-04,
  2.30000000e-04,
  3.40000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
 -1.09752778e-04,
 -2.70745840e-20,
  1.09752778e-04,
  4.48999871e-04,
  5.58752650e-04,
  6.68505428e-04,
 -1.09506112e-04,
 -2.70441593e-20,
  1.09506112e-04,
  4.48001982e-04,
  5.57508094e-04,
  6.67014206e-04,
 -1.09990100e-04,
 -2.71038346e-20,
  1.09990100e-04,
  4.49959952e-04,
  5.59950052e-04,
  6.69940153e-04,
  1.30000000e-04,
  2.40000000e-04,
  3.50000000e-04,
  4.60000000e-04,
  5.70000000e-04,
  6.80000000e-04,
 -1.09987329e-04,
  1.56113122e-24,
  1.09987329e-04,
  4.59947587e-04,
  5.69934916e-04,
  6.79922244e-04,
 -1.09998733e-04,
  1.56123808e-25,
  1.09998733e-04,
  4.59994758e-04,
  5.69993491e-04,
  6.79992224e-04,
 -1.09999873e-04,
  1.56132337e-26,
  1.09999873e-04,
  4.59999476e-04,
  5.69999349e-04,
  6.79999222e-04,
};

pylith::materials::GenMaxwellIsotropic3DTimeDepData::GenMaxwellIsotropic3DTimeDepData(void)
{ // constructor
  dimension = _dimension;
  numLocs = _numLocs;
  numProperties = _numProperties;
  numStateVars = _numStateVars;
  numDBProperties = _numDBProperties;
  numDBStateVars = _numDBStateVars;
  numPropsQuadPt = _numPropsQuadPt;
  numVarsQuadPt = _numVarsQuadPt;
  lengthScale = _lengthScale;
  timeScale = _timeScale;
  pressureScale = _pressureScale;
  densityScale = _densityScale;
  dtStableImplicit = _dtStableImplicit;
  numPropertyValues = const_cast<int*>(_numPropertyValues);
  numStateVarValues = const_cast<int*>(_numStateVarValues);
  dbPropertyValues = const_cast<char**>(_dbPropertyValues);
  dbStateVarValues = const_cast<char**>(_dbStateVarValues);
  dbProperties = const_cast<double*>(_dbProperties);
  dbStateVars = const_cast<double*>(_dbStateVars);
  properties = const_cast<double*>(_properties);
  stateVars = const_cast<double*>(_stateVars);
  propertiesNondim = const_cast<double*>(_propertiesNondim);
  stateVarsNondim = const_cast<double*>(_stateVarsNondim);
  density = const_cast<double*>(_density);
  strain = const_cast<double*>(_strain);
  stress = const_cast<double*>(_stress);
  elasticConsts = const_cast<double*>(_elasticConsts);
  initialStress = const_cast<double*>(_initialStress);
  initialStrain = const_cast<double*>(_initialStrain);
  stateVarsUpdated = const_cast<double*>(_stateVarsUpdated);
} // constructor

pylith::materials::GenMaxwellIsotropic3DTimeDepData::~GenMaxwellIsotropic3DTimeDepData(void)
{}


// End of file
