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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application genmaxwellisotropic3delastic.

#include "GenMaxwellIsotropic3DElasticData.hh"

const int pylith::materials::GenMaxwellIsotropic3DElasticData::_dimension = 3;

const int pylith::materials::GenMaxwellIsotropic3DElasticData::_numLocs = 2;

const int pylith::materials::GenMaxwellIsotropic3DElasticData::_numProperties = 9;

const int pylith::materials::GenMaxwellIsotropic3DElasticData::_numStateVars = 4;

const int pylith::materials::GenMaxwellIsotropic3DElasticData::_numDBProperties = 9;

const int pylith::materials::GenMaxwellIsotropic3DElasticData::_numDBStateVars = 24;

const int pylith::materials::GenMaxwellIsotropic3DElasticData::_numPropsQuadPt = 9;

const int pylith::materials::GenMaxwellIsotropic3DElasticData::_numVarsQuadPt = 24;

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_lengthScale =   1.00000000e+03;

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_timeScale =   1.00000000e+00;

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_pressureScale =   2.25000000e+10;

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_densityScale =   1.00000000e+03;

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_dtStableImplicit =   8.88888889e+06;

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_dtStableExplicit =   1.92450090e-01;

const int pylith::materials::GenMaxwellIsotropic3DElasticData::_numPropertyValues[] = {
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

const int pylith::materials::GenMaxwellIsotropic3DElasticData::_numStateVarValues[] = {
6,
6,
6,
6,
};

const char* pylith::materials::GenMaxwellIsotropic3DElasticData::_dbPropertyValues[] = {
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

const char* pylith::materials::GenMaxwellIsotropic3DElasticData::_dbStateVarValues[] = {
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

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_dbProperties[] = {
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

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_dbStateVars[] = {
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

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_properties[] = {
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

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_stateVars[] = {
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

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_propertiesNondim[] = {
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

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_stateVarsNondim[] = {
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

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_density[] = {
  2.50000000e+03,
  2.00000000e+03,
};

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_strain[] = {
  1.10000000e-04,
  2.20000000e-04,
  3.30000000e-04,
  4.40000000e-04,
  5.50000000e-04,
  6.60000000e-04,
  1.20000000e-04,
  2.30000000e-04,
  3.40000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
};

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_stress[] = {
  1.62660000e+07,
  2.11720000e+07,
  2.60780000e+07,
  1.82940000e+07,
  2.32000000e+07,
  2.81060000e+07,
  1.84236000e+06,
  2.47120000e+06,
  3.10004000e+06,
  2.27736000e+06,
  2.90620000e+06,
  3.53504000e+06,
};

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_elasticConsts[] = {
  6.75000000e+10,
  2.25000000e+10,
  2.25000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.25000000e+10,
  6.75000000e+10,
  2.25000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.25000000e+10,
  2.25000000e+10,
  6.75000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  8.64000000e+09,
  2.88000000e+09,
  2.88000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.88000000e+09,
  8.64000000e+09,
  2.88000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.88000000e+09,
  2.88000000e+09,
  8.64000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
};

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_initialStress[] = {
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

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_initialStrain[] = {
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

const PylithScalar pylith::materials::GenMaxwellIsotropic3DElasticData::_stateVarsUpdated[] = {
  1.10000000e-04,
  2.20000000e-04,
  3.30000000e-04,
  4.40000000e-04,
  5.50000000e-04,
  6.60000000e-04,
 -1.09000000e-04,
  0.00000000e+00,
  1.09000000e-04,
  4.06000000e-04,
  5.15000000e-04,
  6.24000000e-04,
 -1.09000000e-04,
  0.00000000e+00,
  1.09000000e-04,
  4.06000000e-04,
  5.15000000e-04,
  6.24000000e-04,
 -1.09000000e-04,
  0.00000000e+00,
  1.09000000e-04,
  4.06000000e-04,
  5.15000000e-04,
  6.24000000e-04,
  1.20000000e-04,
  2.30000000e-04,
  3.40000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
 -1.09000000e-04,
  2.71050543e-20,
  1.09000000e-04,
  3.86000000e-04,
  4.95000000e-04,
  6.04000000e-04,
 -1.09000000e-04,
  2.71050543e-20,
  1.09000000e-04,
  3.86000000e-04,
  4.95000000e-04,
  6.04000000e-04,
 -1.09000000e-04,
  2.71050543e-20,
  1.09000000e-04,
  3.86000000e-04,
  4.95000000e-04,
  6.04000000e-04,
};

pylith::materials::GenMaxwellIsotropic3DElasticData::GenMaxwellIsotropic3DElasticData(void)
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

pylith::materials::GenMaxwellIsotropic3DElasticData::~GenMaxwellIsotropic3DElasticData(void)
{}


// End of file
