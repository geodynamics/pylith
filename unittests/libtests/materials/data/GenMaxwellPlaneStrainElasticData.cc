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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application genmaxwellplanestrainelastic.

#include "GenMaxwellPlaneStrainElasticData.hh"

const int pylith::materials::GenMaxwellPlaneStrainElasticData::_dimension = 2;

const int pylith::materials::GenMaxwellPlaneStrainElasticData::_numLocs = 2;

const int pylith::materials::GenMaxwellPlaneStrainElasticData::_numProperties = 9;

const int pylith::materials::GenMaxwellPlaneStrainElasticData::_numStateVars = 5;

const int pylith::materials::GenMaxwellPlaneStrainElasticData::_numDBProperties = 9;

const int pylith::materials::GenMaxwellPlaneStrainElasticData::_numDBStateVars = 16;

const int pylith::materials::GenMaxwellPlaneStrainElasticData::_numPropsQuadPt = 9;

const int pylith::materials::GenMaxwellPlaneStrainElasticData::_numVarsQuadPt = 16;

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_lengthScale =   1.00000000e+03;

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_timeScale =   1.00000000e+00;

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_pressureScale =   2.25000000e+10;

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_densityScale =   1.00000000e+03;

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_dtStableImplicit =   8.88888889e+06;

const int pylith::materials::GenMaxwellPlaneStrainElasticData::_numPropertyValues[] = {
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

const int pylith::materials::GenMaxwellPlaneStrainElasticData::_numStateVarValues[] = {
1,
3,
4,
4,
4,
};

const char* pylith::materials::GenMaxwellPlaneStrainElasticData::_dbPropertyValues[] = {
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

const char* pylith::materials::GenMaxwellPlaneStrainElasticData::_dbStateVarValues[] = {
"stress-zz-initial",
"total-strain-xx",
"total-strain-yy",
"total-strain-xy",
"viscous-strain-1-xx",
"viscous-strain-1-yy",
"viscous-strain-1-zz",
"viscous-strain-1-xy",
"viscous-strain-2-xx",
"viscous-strain-2-yy",
"viscous-strain-2-zz",
"viscous-strain-2-xy",
"viscous-strain-3-xx",
"viscous-strain-3-yy",
"viscous-strain-3-zz",
"viscous-strain-3-xy",
};

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_dbProperties[] = {
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

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_dbStateVars[] = {
  1.50000000e+04,
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
  4.50000000e+04,
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

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_properties[] = {
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

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_stateVars[] = {
  1.50000000e+04,
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
  4.50000000e+04,
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

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_propertiesNondim[] = {
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

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_stateVarsNondim[] = {
  6.66666667e-07,
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
  2.00000000e-06,
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

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_density[] = {
  2.50000000e+03,
  2.00000000e+03,
};

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_strain[] = {
  1.10000000e-04,
  1.20000000e-04,
  1.40000000e-04,
  4.10000000e-04,
  4.20000000e-04,
  4.40000000e-04,
};

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_stress[] = {
  7.33350000e+06,
  7.73950000e+06,
  4.79400000e+06,
  4.09740000e+06,
  4.15024000e+06,
  2.21976000e+06,
};

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_elasticConsts[] = {
  6.75000000e+10,
  2.25000000e+10,
  0.00000000e+00,
  2.25000000e+10,
  6.75000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  8.64000000e+09,
  2.88000000e+09,
  0.00000000e+00,
  2.88000000e+09,
  8.64000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
};

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_initialStress[] = {
  2.10000000e+04,
  2.20000000e+04,
  2.40000000e+04,
  5.10000000e+04,
  5.20000000e+04,
  5.40000000e+04,
};

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_initialStrain[] = {
  3.10000000e-05,
  3.20000000e-05,
  3.40000000e-05,
  6.10000000e-05,
  6.20000000e-05,
  6.40000000e-05,
};

const double pylith::materials::GenMaxwellPlaneStrainElasticData::_stateVarsUpdated[] = {
  1.50000000e+04,
  1.10000000e-04,
  1.20000000e-04,
  1.40000000e-04,
  2.33333333e-05,
  3.23333333e-05,
 -5.56666667e-05,
  1.06000000e-04,
  2.33333333e-05,
  3.23333333e-05,
 -5.56666667e-05,
  1.06000000e-04,
  2.33333333e-05,
  3.23333333e-05,
 -5.56666667e-05,
  1.06000000e-04,
  4.50000000e+04,
  4.10000000e-04,
  4.20000000e-04,
  4.40000000e-04,
  1.13333333e-04,
  1.22333333e-04,
 -2.35666667e-04,
  3.76000000e-04,
  1.13333333e-04,
  1.22333333e-04,
 -2.35666667e-04,
  3.76000000e-04,
  1.13333333e-04,
  1.22333333e-04,
 -2.35666667e-04,
  3.76000000e-04,
};

pylith::materials::GenMaxwellPlaneStrainElasticData::GenMaxwellPlaneStrainElasticData(void)
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

pylith::materials::GenMaxwellPlaneStrainElasticData::~GenMaxwellPlaneStrainElasticData(void)
{}


// End of file
