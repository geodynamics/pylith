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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application druckerprager3dtimedep.

#include "DruckerPrager3DTimeDepData.hh"

const int pylith::materials::DruckerPrager3DTimeDepData::_dimension = 3;

const int pylith::materials::DruckerPrager3DTimeDepData::_numLocs = 2;

const int pylith::materials::DruckerPrager3DTimeDepData::_numProperties = 6;

const int pylith::materials::DruckerPrager3DTimeDepData::_numStateVars = 1;

const int pylith::materials::DruckerPrager3DTimeDepData::_numDBProperties = 6;

const int pylith::materials::DruckerPrager3DTimeDepData::_numDBStateVars = 6;

const int pylith::materials::DruckerPrager3DTimeDepData::_numPropsQuadPt = 6;

const int pylith::materials::DruckerPrager3DTimeDepData::_numVarsQuadPt = 6;

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_lengthScale =   1.00000000e+03;

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_timeScale =   1.00000000e+00;

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_pressureScale =   2.25000000e+10;

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_densityScale =   1.00000000e+03;

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_dtStableImplicit =   1.00000000e+99;

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_dtStableExplicit =   1.92450090e-01;

const int pylith::materials::DruckerPrager3DTimeDepData::_numPropertyValues[] = {
1,
1,
1,
1,
1,
1,
};

const int pylith::materials::DruckerPrager3DTimeDepData::_numStateVarValues[] = {
6,
};

const char* pylith::materials::DruckerPrager3DTimeDepData::_dbPropertyValues[] = {
"density",
"vs",
"vp",
"friction-angle",
"cohesion",
"dilatation-angle",
};

const char* pylith::materials::DruckerPrager3DTimeDepData::_dbStateVarValues[] = {
"plastic-strain-xx",
"plastic-strain-yy",
"plastic-strain-zz",
"plastic-strain-xy",
"plastic-strain-yz",
"plastic-strain-xz",
};

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_dbProperties[] = {
  2.50000000e+03,
  3.00000000e+03,
  5.19615242e+03,
  5.23598776e-01,
  3.00000000e+05,
  3.49065850e-01,
  2.00000000e+03,
  1.20000000e+03,
  2.07846097e+03,
  4.36332313e-01,
  1.00000000e+04,
  4.36332313e-01,
};

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_dbStateVars[] = {
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

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_properties[] = {
  2.50000000e+03,
  2.25000000e+10,
  2.25000000e+10,
  2.30940108e-01,
  3.60000000e+05,
  1.48583084e-01,
  2.00000000e+03,
  2.88000000e+09,
  2.88000000e+09,
  1.89338478e-01,
  1.21811303e+04,
  1.89338478e-01,
};

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_stateVars[] = {
  4.10000000e-05,
  4.20000000e-05,
  4.30000000e-05,
  4.40000000e-05,
  4.50000000e-05,
  4.60000000e-05,
  1.10000000e-05,
  1.20000000e-05,
  1.30000000e-05,
  1.40000000e-05,
  1.50000000e-05,
  1.60000000e-05,
};

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_propertiesNondim[] = {
  2.50000000e+00,
  1.00000000e+00,
  1.00000000e+00,
  2.30940108e-01,
  1.60000000e-05,
  1.48583084e-01,
  2.00000000e+00,
  1.28000000e-01,
  1.28000000e-01,
  1.89338478e-01,
  5.41383567e-07,
  1.89338478e-01,
};

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_stateVarsNondim[] = {
  4.10000000e-05,
  4.20000000e-05,
  4.30000000e-05,
  4.40000000e-05,
  4.50000000e-05,
  4.60000000e-05,
  1.10000000e-05,
  1.20000000e-05,
  1.30000000e-05,
  1.40000000e-05,
  1.50000000e-05,
  1.60000000e-05,
};

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_density[] = {
  2.50000000e+03,
  2.00000000e+03,
};

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_strain[] = {
 -2.10000000e-04,
  1.20000000e-04,
  1.30000000e-04,
  1.10000000e-05,
  1.10000000e-05,
  1.10000000e-05,
  4.10000000e-04,
  4.20000000e-04,
  4.30000000e-04,
  1.40000000e-04,
  1.50000000e-04,
  1.60000000e-04,
};

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_stress[] = {
 -1.63946744e+07,
 -6.17719414e+06,
 -5.86690637e+06,
 -2.02684646e+06,
 -2.02615846e+06,
 -2.02547046e+06,
  4.21392632e+06,
  4.21392632e+06,
  4.21392632e+06,
  1.45519152e-11,
  7.27595761e-12,
  0.00000000e+00,
};

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_elasticConsts[] = {
  6.57004938e+10,
  2.99833221e+10,
  2.98388535e+10,
  1.88737456e+09,
  1.88673474e+09,
  1.88609213e+09,
  2.09815549e+10,
  3.02171782e+10,
 -1.40253222e+09,
  8.61894060e+09,
  8.61601438e+09,
  8.61308910e+09,
  2.05637179e+10,
 -1.67589867e+09,
  2.86086951e+10,
  8.82336730e+09,
  8.82037124e+09,
  8.81737750e+09,
  2.72937212e+09,
  6.09515538e+09,
  6.19736873e+09,
  2.96246322e+10,
 -1.33489072e+09,
 -1.33443763e+09,
  2.72844546e+09,
  6.09308621e+09,
  6.19526533e+09,
 -1.33489084e+09,
  2.96255386e+10,
 -1.33398466e+09,
  2.72751914e+09,
  6.09101739e+09,
  6.19316124e+09,
 -1.33443763e+09,
 -1.33398443e+09,
  2.96264446e+10,
  4.98723239e+09,
  4.80000023e+09,
  4.61276807e+09,
 -2.75144819e+09,
 -3.12591158e+09,
 -3.50037543e+09,
  4.98723239e+09,
  4.80000023e+09,
  4.61276807e+09,
 -2.75144819e+09,
 -3.12591158e+09,
 -3.50037543e+09,
  4.98723239e+09,
  4.80000023e+09,
  4.61276807e+09,
 -2.75144819e+09,
 -3.12591158e+09,
 -3.50037543e+09,
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

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_initialStress[] = {
  2.10000000e+04,
  2.20000000e+04,
  2.30000000e+04,
  2.40000000e+04,
  2.50000000e+04,
  2.60000000e+04,
  5.60000000e+04,
  5.50000000e+04,
  5.40000000e+04,
  5.30000000e+04,
  5.20000000e+04,
  5.10000000e+04,
};

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_initialStrain[] = {
  3.60000000e-05,
  3.50000000e-05,
  3.40000000e-05,
  3.30000000e-05,
  3.20000000e-05,
  3.10000000e-05,
  6.60000000e-05,
  6.50000000e-05,
  6.40000000e-05,
  6.30000000e-05,
  6.20000000e-05,
  6.10000000e-05,
};

const PylithScalar pylith::materials::DruckerPrager3DTimeDepData::_stateVarsUpdated[] = {
 -7.89512481e-06,
  9.60719814e-05,
  1.00198920e-04,
  2.35743657e-05,
  2.45812990e-05,
  2.55882324e-05,
  5.53592834e-05,
  6.61856723e-05,
  7.70120612e-05,
  8.62013889e-05,
  9.70277778e-05,
  1.07854167e-04,
};

pylith::materials::DruckerPrager3DTimeDepData::DruckerPrager3DTimeDepData(void)
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
  dtStableExplicit = _dtStableExplicit;
  numPropertyValues = const_cast<int*>(_numPropertyValues);
  numStateVarValues = const_cast<int*>(_numStateVarValues);
  dbPropertyValues = const_cast<char**>(_dbPropertyValues);
  dbStateVarValues = const_cast<char**>(_dbStateVarValues);
  dbProperties = const_cast<PylithScalar*>(_dbProperties);
  dbStateVars = const_cast<PylithScalar*>(_dbStateVars);
  properties = const_cast<PylithScalar*>(_properties);
  stateVars = const_cast<PylithScalar*>(_stateVars);
  propertiesNondim = const_cast<PylithScalar*>(_propertiesNondim);
  stateVarsNondim = const_cast<PylithScalar*>(_stateVarsNondim);
  density = const_cast<PylithScalar*>(_density);
  strain = const_cast<PylithScalar*>(_strain);
  stress = const_cast<PylithScalar*>(_stress);
  elasticConsts = const_cast<PylithScalar*>(_elasticConsts);
  initialStress = const_cast<PylithScalar*>(_initialStress);
  initialStrain = const_cast<PylithScalar*>(_initialStrain);
  stateVarsUpdated = const_cast<PylithScalar*>(_stateVarsUpdated);
} // constructor

pylith::materials::DruckerPrager3DTimeDepData::~DruckerPrager3DTimeDepData(void)
{}


// End of file
