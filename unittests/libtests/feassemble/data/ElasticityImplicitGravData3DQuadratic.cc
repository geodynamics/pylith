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
// This file was generated from python application elasticityapp.

#include "ElasticityImplicitGravData3DQuadratic.hh"

const int pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_spaceDim = 3;

const int pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_cellDim = 3;

const int pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_numVertices = 10;

const int pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_numCells = 1;

const int pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_numBasis = 10;

const int pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_numQuadPts = 4;

const char* pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_matType = "ElasticIsotropic3D";

const char* pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_matDBFilename = "data/elasticisotropic3d.spatialdb";

const int pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_matId = 0;

const char* pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_matLabel = "elastic isotropic 3-D";

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_dt =   1.00000000e-02;

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_gravityVec[] = {
  0.00000000e+00,  0.00000000e+00, -1.00000000e+08,
};

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_vertices[] = {
 -5.00000000e-01, -2.00000000e+00, -1.00000000e+00,
  2.00000000e+00, -2.00000000e+00, -5.00000000e-01,
  1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
  2.00000000e-01,  5.00000000e-01,  2.00000000e+00,
  1.50000000e+00, -5.00000000e-01, -2.50000000e-01,
  2.50000000e-01, -5.00000000e-01, -5.00000000e-01,
  7.50000000e-01, -2.00000000e+00, -7.50000000e-01,
 -1.50000000e-01, -7.50000000e-01,  5.00000000e-01,
  1.10000000e+00, -7.50000000e-01,  7.50000000e-01,
  6.00000000e-01,  7.50000000e-01,  1.00000000e+00,
};

const int pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_cells[] = {
0,1,2,3,4,5,6,7,8,9,
};

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,  1.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  0.00000000e+00, -1.00000000e+00,
  0.00000000e+00, -1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -1.00000000e+00,  0.00000000e+00,
 -1.00000000e+00,  0.00000000e+00,  0.00000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_quadPts[] = {
 -8.00000000e-01, -8.00000000e-01, -8.00000000e-01,
  5.00000000e-01, -8.00000000e-01, -8.00000000e-01,
 -8.00000000e-01,  5.00000000e-01, -8.00000000e-01,
 -8.00000000e-01, -8.00000000e-01,  5.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_quadWts[] = {
  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_basis[] = {
  2.80000000e-01, -8.00000000e-02, -8.00000000e-02,
 -8.00000000e-02,  4.00000000e-02,  2.80000000e-01,
  2.80000000e-01,  2.80000000e-01,  4.00000000e-02,
  4.00000000e-02, -4.50000000e-02,  3.75000000e-01,
 -8.00000000e-02, -8.00000000e-02,  3.00000000e-01,
  2.00000000e-02,  1.50000000e-01,  2.00000000e-02,
  3.00000000e-01,  4.00000000e-02, -4.50000000e-02,
 -8.00000000e-02,  3.75000000e-01, -8.00000000e-02,
  3.00000000e-01,  1.50000000e-01,  2.00000000e-02,
  2.00000000e-02,  4.00000000e-02,  3.00000000e-01,
 -4.50000000e-02, -8.00000000e-02, -8.00000000e-02,
  3.75000000e-01,  4.00000000e-02,  2.00000000e-02,
  2.00000000e-02,  1.50000000e-01,  3.00000000e-01,
  3.00000000e-01,};

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_basisDerivRef[] = {
 -9.00000000e-01, -9.00000000e-01, -9.00000000e-01,
 -3.00000000e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -3.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -3.00000000e-01,
  2.00000000e-01,  2.00000000e-01,  0.00000000e+00,
 -2.00000000e-01,  1.20000000e+00, -2.00000000e-01,
  1.20000000e+00, -2.00000000e-01, -2.00000000e-01,
 -2.00000000e-01, -2.00000000e-01,  1.20000000e+00,
  2.00000000e-01,  0.00000000e+00,  2.00000000e-01,
  0.00000000e+00,  2.00000000e-01,  2.00000000e-01,
  4.00000000e-01,  4.00000000e-01,  4.00000000e-01,
  1.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -3.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -3.00000000e-01,
  2.00000000e-01,  1.50000000e+00,  0.00000000e+00,
 -2.00000000e-01, -1.00000000e-01, -2.00000000e-01,
 -1.40000000e+00, -1.50000000e+00, -1.50000000e+00,
 -2.00000000e-01, -2.00000000e-01, -1.00000000e-01,
  2.00000000e-01,  0.00000000e+00,  1.50000000e+00,
  0.00000000e+00,  2.00000000e-01,  2.00000000e-01,
  4.00000000e-01,  4.00000000e-01,  4.00000000e-01,
 -3.00000000e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -3.00000000e-01,
  1.50000000e+00,  2.00000000e-01,  0.00000000e+00,
 -1.50000000e+00, -1.40000000e+00, -1.50000000e+00,
 -1.00000000e-01, -2.00000000e-01, -2.00000000e-01,
 -2.00000000e-01, -2.00000000e-01, -1.00000000e-01,
  2.00000000e-01,  0.00000000e+00,  2.00000000e-01,
  0.00000000e+00,  2.00000000e-01,  1.50000000e+00,
  4.00000000e-01,  4.00000000e-01,  4.00000000e-01,
 -3.00000000e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -3.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.00000000e+00,
  2.00000000e-01,  2.00000000e-01,  0.00000000e+00,
 -2.00000000e-01, -1.00000000e-01, -2.00000000e-01,
 -1.00000000e-01, -2.00000000e-01, -2.00000000e-01,
 -1.50000000e+00, -1.50000000e+00, -1.40000000e+00,
  1.50000000e+00,  0.00000000e+00,  2.00000000e-01,
  0.00000000e+00,  1.50000000e+00,  2.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_fieldTIncr[] = {
  3.00000000e-01, -4.00000000e-01, -4.00000000e-01,
 -6.00000000e-01,  8.00000000e-01,  2.00000000e-01,
  5.00000000e-01,  5.00000000e-01,  7.00000000e-01,
 -7.00000000e-01, -5.00000000e-01, -7.00000000e-01,
 -6.00000000e-01, -3.00000000e-01,  8.00000000e-01,
 -4.00000000e-01, -8.00000000e-01, -5.00000000e-01,
  7.00000000e-01,  8.00000000e-01, -5.00000000e-01,
 -5.00000000e-01, -5.00000000e-01, -7.00000000e-01,
 -3.00000000e-01, -9.00000000e-01,  8.00000000e-01,
 -1.00000000e-01,  5.00000000e-01, -9.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_fieldT[] = {
  1.00000000e-01, -2.00000000e-01, -6.00000000e-01,
 -3.00000000e-01,  4.00000000e-01,  9.00000000e-01,
  6.00000000e-01,  8.00000000e-01,  5.00000000e-01,
 -8.00000000e-01, -6.00000000e-01, -8.00000000e-01,
 -0.00000000e+00, -2.00000000e-01,  6.00000000e-01,
 -4.00000000e-01, -7.00000000e-01, -2.00000000e-01,
  7.00000000e-01,  6.00000000e-01, -1.00000000e-01,
 -4.00000000e-01, -3.00000000e-01, -3.00000000e-01,
 -7.00000000e-01, -6.00000000e-01,  1.00000000e-01,
 -9.00000000e-01,  3.00000000e-01, -8.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_fieldTmdt[] = {
  2.00000000e-01, -3.00000000e-01, -1.00000000e-01,
 -4.00000000e-01,  2.00000000e-01,  3.00000000e-01,
 -5.00000000e-01,  2.00000000e-01,  2.00000000e-01,
 -3.00000000e-01, -8.00000000e-01, -3.00000000e-01,
 -5.00000000e-01, -9.00000000e-01,  4.00000000e-01,
 -3.00000000e-01, -6.00000000e-01, -8.00000000e-01,
  9.00000000e-01,  5.00000000e-01, -2.00000000e-01,
 -7.00000000e-01, -2.00000000e-01, -9.00000000e-01,
 -5.00000000e-01, -8.00000000e-01,  4.00000000e-01,
 -4.00000000e-01,  5.00000000e-01, -7.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_valsResidual[] = {
  2.17091508e+10, -1.76856611e+10, -3.16592972e+10,
  1.31632837e+11, -4.68038741e+10, -1.57191521e+10,
 -5.63686735e+10, -5.92797533e+10, -1.37382956e+11,
 -5.78258053e+09,  1.17383971e+11, -1.67395324e+11,
  5.81040000e+10,  1.11498152e+11, -1.28438506e+11,
  1.14913086e+11,  1.99398526e+11, -8.74985277e+10,
 -1.78723310e+11, -1.16643455e+11, -8.18170460e+10,
  6.28557950e+10, -9.05944297e+10,  8.36876150e+09,
 -1.94927791e+11,  1.49434559e+11, -2.43345729e+11,
  4.65874861e+10, -2.46708036e+11,  1.73429444e+11,
};

const PylithScalar pylith::feassemble::ElasticityImplicitGravData3DQuadratic::_valsJacobian[] = {
  4.84108858e+10,  7.55490483e+09,  1.36932650e+10,
  1.51828404e+10, -1.68096633e+09,  3.67199854e+09,
  2.10939239e+09,  5.87866032e+09, -4.46211567e+09,
 -1.15527086e+09, -1.67939239e+09,  5.35453880e+09,
  2.33244070e+10,  5.66200586e+09, -1.06573939e+09,
 -2.73584773e+10, -1.89820644e+10,  5.67313324e+09,
 -6.20183163e+10,  1.05973646e+09, -1.58917277e+10,
 -1.87033236e+10,  1.05556369e+09, -2.03524158e+10,
  1.89209078e+10, -4.53257687e+09,  1.21753294e+10,
  1.28695461e+09,  5.66412884e+09,  1.20373353e+09,
  7.55490483e+09,  2.08335944e+10,  3.50547584e+09,
 -1.68096633e+09,  3.42111274e+09, -1.05516837e+09,
  5.87866032e+09,  3.02290630e+09,  1.35169839e+09,
 -1.67939239e+09,  5.00512445e+08,  8.71961933e+08,
  5.66200586e+09,  8.69193265e+09,  3.99970717e+08,
 -1.89820644e+10, -1.73812592e+10, -5.15967789e+09,
  1.05973646e+09, -1.84369693e+10,  1.22131772e+09,
  1.05556369e+09, -1.06939824e+10, -3.88781845e+09,
 -4.53257687e+09,  5.28963397e+09, -2.47115666e+08,
  5.66412884e+09,  4.75251830e+09,  2.99935578e+09,
  1.36932650e+10,  3.50547584e+09,  2.52532138e+10,
  3.67199854e+09, -1.05516837e+09,  6.15724744e+09,
 -4.46211567e+09,  1.35169839e+09, -4.42245242e+09,
  5.35453880e+09,  8.71961933e+08,  6.68294290e+09,
 -1.06573939e+09,  3.99970717e+08,  2.33995608e+09,
  5.67313324e+09, -5.15967789e+09,  3.70483163e+08,
 -1.58917277e+10,  1.22131772e+09, -2.76780234e+10,
 -2.03524158e+10, -3.88781845e+09, -2.90717277e+10,
  1.21753294e+10, -2.47115666e+08,  1.73193265e+10,
  1.20373353e+09,  2.99935578e+09,  3.04903367e+09,
  1.51828404e+10, -1.68096633e+09,  3.67199854e+09,
  5.35194363e+10, -1.83713031e+10,  7.97701318e+09,
  1.61204978e+09, -5.13162518e+09,  4.74267936e+09,
  1.25754026e+09,  6.97950220e+08, -5.69121523e+09,
 -2.06587189e+10,  2.14802343e+10, -1.72295754e+10,
  1.45228258e+10, -6.08931186e+09,  7.45062225e+09,
 -6.26214788e+10,  1.06571010e+10, -1.36590044e+10,
  1.42020791e+10, -8.14934114e+08, -1.98956808e+09,
 -1.95445461e+10,  3.15871157e+09,  1.55626647e+10,
  2.52797218e+09, -3.90585652e+09, -8.35614934e+08,
 -1.68096633e+09,  3.42111274e+09, -1.05516837e+09,
 -1.83713031e+10,  3.28349854e+10, -4.66348463e+09,
 -5.13162518e+09,  1.02936896e+10, -4.37187408e+09,
  6.97950220e+08, -2.77642753e+09,  3.86024890e+09,
  2.14802343e+10, -4.19176940e+10,  1.50988287e+10,
 -6.08931186e+09,  1.22570864e+10, -4.86344070e+09,
  1.06571010e+10, -2.03266618e+10,  4.63449488e+09,
 -8.14934114e+08,  4.31742313e+08,  2.58467057e+09,
  3.15871157e+09, -8.40183016e+08, -1.07735578e+10,
 -3.90585652e+09,  6.62234993e+09, -4.50717423e+08,
  3.67199854e+09, -1.05516837e+09,  6.15724744e+09,
  7.97701318e+09, -4.66348463e+09,  2.41197731e+10,
  4.74267936e+09, -4.37187408e+09,  6.14718887e+09,
 -5.69121523e+09,  3.86024890e+09, -4.18462665e+09,
 -1.72295754e+10,  1.50988287e+10, -2.63467570e+10,
  7.45062225e+09, -4.86344070e+09,  1.08598389e+10,
 -1.36590044e+10,  4.63449488e+09, -2.61182577e+10,
 -1.98956808e+09,  2.58467057e+09,  1.51200586e+09,
  1.55626647e+10, -1.07735578e+10,  6.12466325e+09,
 -8.35614934e+08, -4.50717423e+08,  1.72892387e+09,
  2.10939239e+09,  5.87866032e+09, -4.46211567e+09,
  1.61204978e+09, -5.13162518e+09,  4.74267936e+09,
  3.52820132e+10,  6.64751098e+09, -5.81076135e+09,
  7.99568082e+09,  1.58806735e+09, -2.30600293e+09,
 -1.53454539e+10,  1.41912884e+10, -1.32126647e+10,
 -1.70321669e+10, -2.00351391e+10,  1.54363104e+10,
  3.27357980e+09,  4.15483163e+08,  4.51500732e+08,
  9.04924597e+09,  6.49520498e+09, -5.92587848e+09,
  8.46395315e+09, -3.12170571e+09,  2.14659590e+09,
 -3.54082943e+10, -6.92774524e+09,  8.94033675e+09,
  5.87866032e+09,  3.02290630e+09,  1.35169839e+09,
 -5.13162518e+09,  1.02936896e+10, -4.37187408e+09,
  6.64751098e+09,  7.21436237e+10, -3.32375549e+10,
  1.58806735e+09,  1.06122255e+10, -7.94033675e+09,
  1.41912884e+10, -5.33698170e+10,  2.34235578e+10,
 -2.00351391e+10, -3.08659004e+10,  5.79569546e+09,
  4.15483163e+08,  1.19144436e+10, -2.79241581e+09,
  6.49520498e+09,  1.22026428e+10, -6.02102489e+09,
 -3.12170571e+09,  1.84171157e+10, -1.08464714e+10,
 -6.92774524e+09, -5.43709297e+10,  3.46387262e+10,
 -4.46211567e+09,  1.35169839e+09, -4.42245242e+09,
  4.74267936e+09, -4.37187408e+09,  6.14718887e+09,
 -5.81076135e+09, -3.32375549e+10,  6.31736676e+10,
 -2.30600293e+09, -7.94033675e+09,  1.90644949e+10,
 -1.32126647e+10,  2.34235578e+10, -3.77246633e+10,
  1.54363104e+10,  5.79569546e+09, -5.32645681e+09,
  4.51500732e+08, -2.79241581e+09,  1.75637628e+09,
 -5.92587848e+09, -6.02102489e+09,  1.34434627e+10,
  2.14659590e+09, -1.08464714e+10,  2.22102928e+10,
  8.94033675e+09,  3.46387262e+10, -7.83219107e+10,
 -1.15527086e+09, -1.67939239e+09,  5.35453880e+09,
  1.25754026e+09,  6.97950220e+08, -5.69121523e+09,
  7.99568082e+09,  1.58806735e+09, -2.30600293e+09,
  2.45678990e+10,  1.95241581e+09, -8.36749634e+09,
  8.15164714e+09,  2.01387262e+09, -7.04516837e+09,
  6.24011713e+09, -8.23572474e+06,  2.52101025e+09,
  1.43704246e+08, -8.13579795e+08, -5.41800878e+08,
 -3.80461201e+09,  4.57115666e+09, -1.39335286e+10,
 -1.11098463e+10, -2.76237189e+09,  2.03244510e+10,
 -3.22868594e+10, -5.55988287e+09,  9.68521230e+09,
 -1.67939239e+09,  5.00512445e+08,  8.71961933e+08,
  6.97950220e+08, -2.77642753e+09,  3.86024890e+09,
  1.58806735e+09,  1.06122255e+10, -7.94033675e+09,
  1.95241581e+09,  2.51722182e+10, -9.76207906e+09,
  2.01387262e+09,  6.90296486e+09, -3.59436310e+09,
 -8.23572474e+06,  1.00322108e+10, -6.43382138e+09,
 -8.13579795e+08, -2.08133236e+09,  4.24289898e+09,
  4.57115666e+09, -9.06830161e+09,  2.44216691e+08,
 -2.76237189e+09,  1.39227672e+09, -9.28814056e+09,
 -5.55988287e+09, -4.06863470e+10,  2.77994143e+10,
  5.35453880e+09,  8.71961933e+08,  6.68294290e+09,
 -5.69121523e+09,  3.86024890e+09, -4.18462665e+09,
 -2.30600293e+09, -7.94033675e+09,  1.90644949e+10,
 -8.36749634e+09, -9.76207906e+09,  6.47318814e+10,
 -7.04516837e+09, -3.59436310e+09,  1.31084553e+10,
  2.52101025e+09, -6.43382138e+09,  2.29992679e+10,
 -5.41800878e+08,  4.24289898e+09,  1.96434846e+09,
 -1.39335286e+10,  2.44216691e+08, -3.98836750e+10,
  2.03244510e+10, -9.28814056e+09, -5.70721083e+09,
  9.68521230e+09,  2.77994143e+10, -7.87758785e+10,
  2.33244070e+10,  5.66200586e+09, -1.06573939e+09,
 -2.06587189e+10,  2.14802343e+10, -1.72295754e+10,
 -1.53454539e+10,  1.41912884e+10, -1.32126647e+10,
  8.15164714e+09,  2.01387262e+09, -7.04516837e+09,
  1.60496120e+11, -5.26156662e+09, -1.13147145e+10,
 -7.57647877e+10,  8.17101025e+09, -2.06400439e+10,
 -1.76338141e+10, -3.07990117e+10,  1.64314788e+10,
 -3.77006515e+10, -9.19637628e+09,  1.01987555e+10,
 -2.11084114e+10, -1.16195095e+10,  1.88193265e+10,
 -3.76033675e+09,  5.35805271e+09,  2.50583455e+10,
  5.66200586e+09,  8.69193265e+09,  3.99970717e+08,
  2.14802343e+10, -4.19176940e+10,  1.50988287e+10,
  1.41912884e+10, -5.33698170e+10,  2.34235578e+10,
  2.01387262e+09,  6.90296486e+09, -3.59436310e+09,
 -5.26156662e+09,  1.62573016e+11, -5.65721669e+10,
  8.17101025e+09, -8.45661786e+09, -3.73005124e+09,
 -3.07990117e+10, -2.14388653e+10, -2.61494143e+09,
 -9.19637628e+09, -1.89591728e+10,  4.09688141e+09,
 -1.16195095e+10, -5.44731406e+10,  4.31375476e+10,
  5.35805271e+09,  2.04473939e+10, -1.96452635e+10,
 -1.06573939e+09,  3.99970717e+08,  2.33995608e+09,
 -1.72295754e+10,  1.50988287e+10, -2.63467570e+10,
 -1.32126647e+10,  2.34235578e+10, -3.77246633e+10,
 -7.04516837e+09, -3.59436310e+09,  1.31084553e+10,
 -1.13147145e+10, -5.65721669e+10,  1.42996750e+11,
 -2.06400439e+10, -3.73005124e+09, -1.04425769e+10,
  1.64314788e+10, -2.61494143e+09,  1.66850878e+10,
  1.01987555e+10,  4.09688141e+09, -1.93946779e+10,
  1.88193265e+10,  4.31375476e+10, -9.78411786e+10,
  2.50583455e+10, -1.96452635e+10,  1.66196047e+10,
 -2.73584773e+10, -1.89820644e+10,  5.67313324e+09,
  1.45228258e+10, -6.08931186e+09,  7.45062225e+09,
 -1.70321669e+10, -2.00351391e+10,  1.54363104e+10,
  6.24011713e+09, -8.23572474e+06,  2.52101025e+09,
 -7.57647877e+10,  8.17101025e+09, -2.06400439e+10,
  1.39614114e+11,  1.14450952e+10,  2.58067350e+10,
 -1.01994143e+09,  2.37083821e+10, -1.43894949e+10,
 -1.34082577e+10, -8.13188141e+09,  1.71687775e+10,
 -3.29346413e+10,  9.27946559e+09, -1.70012811e+10,
  7.14121523e+09,  6.42679356e+08, -2.20257687e+10,
 -1.89820644e+10, -1.73812592e+10, -5.15967789e+09,
 -6.08931186e+09,  1.22570864e+10, -4.86344070e+09,
 -2.00351391e+10, -3.08659004e+10,  5.79569546e+09,
 -8.23572474e+06,  1.00322108e+10, -6.43382138e+09,
  8.17101025e+09, -8.45661786e+09, -3.73005124e+09,
  1.14450952e+10,  1.44711406e+11, -4.85954758e+10,
  2.37083821e+10, -4.48338946e+10,  2.22380893e+10,
 -8.13188141e+09, -4.31508638e+10,  3.70394070e+10,
  9.27946559e+09, -2.95963543e+10,  1.40676720e+10,
  6.42679356e+08,  7.28418741e+09, -1.03583968e+10,
  5.67313324e+09, -5.15967789e+09,  3.70483163e+08,
  7.45062225e+09, -4.86344070e+09,  1.08598389e+10,
  1.54363104e+10,  5.79569546e+09, -5.32645681e+09,
  2.52101025e+09, -6.43382138e+09,  2.29992679e+10,
 -2.06400439e+10, -3.73005124e+09, -1.04425769e+10,
  2.58067350e+10, -4.85954758e+10,  1.55051786e+11,
 -1.43894949e+10,  2.22380893e+10, -1.77103660e+10,
  1.71687775e+10,  3.70394070e+10, -7.94883895e+10,
 -1.70012811e+10,  1.40676720e+10, -4.85184919e+10,
 -2.20257687e+10, -1.03583968e+10, -2.77950952e+10,
 -6.20183163e+10,  1.05973646e+09, -1.58917277e+10,
 -6.26214788e+10,  1.06571010e+10, -1.36590044e+10,
  3.27357980e+09,  4.15483163e+08,  4.51500732e+08,
  1.43704246e+08, -8.13579795e+08, -5.41800878e+08,
 -1.76338141e+10, -3.07990117e+10,  1.64314788e+10,
 -1.01994143e+09,  2.37083821e+10, -1.43894949e+10,
  1.39611698e+11, -1.12319180e+10,  2.96296486e+10,
  4.22811127e+09,  4.34842606e+09,  2.43174597e+10,
  7.56508053e+08,  4.33459004e+09, -2.59982430e+10,
 -4.72005124e+09, -1.67920937e+09, -3.49816984e+08,
  1.05973646e+09, -1.84369693e+10,  1.22131772e+09,
  1.06571010e+10, -2.03266618e+10,  4.63449488e+09,
  4.15483163e+08,  1.19144436e+10, -2.79241581e+09,
 -8.13579795e+08, -2.08133236e+09,  4.24289898e+09,
 -3.07990117e+10, -2.14388653e+10, -2.61494143e+09,
  2.37083821e+10, -4.48338946e+10,  2.22380893e+10,
 -1.12319180e+10,  8.22470571e+10,  1.43959004e+09,
  4.34842606e+09,  1.98506003e+10, -1.62421303e+10,
  4.33459004e+09,  7.08971449e+09, -9.18295022e+09,
 -1.67920937e+09, -1.39840922e+10, -2.94395315e+09,
 -1.58917277e+10,  1.22131772e+09, -2.76780234e+10,
 -1.36590044e+10,  4.63449488e+09, -2.61182577e+10,
  4.51500732e+08, -2.79241581e+09,  1.75637628e+09,
 -5.41800878e+08,  4.24289898e+09,  1.96434846e+09,
  1.64314788e+10, -2.61494143e+09,  1.66850878e+10,
 -1.43894949e+10,  2.22380893e+10, -1.77103660e+10,
  2.96296486e+10,  1.43959004e+09,  7.07293851e+10,
  2.43174597e+10, -1.62421303e+10,  1.91043045e+10,
 -2.59982430e+10, -9.18295022e+09, -3.29619253e+10,
 -3.49816984e+08, -2.94395315e+09, -5.77092972e+09,
 -1.87033236e+10,  1.05556369e+09, -2.03524158e+10,
  1.42020791e+10, -8.14934114e+08, -1.98956808e+09,
  9.04924597e+09,  6.49520498e+09, -5.92587848e+09,
 -3.80461201e+09,  4.57115666e+09, -1.39335286e+10,
 -3.77006515e+10, -9.19637628e+09,  1.01987555e+10,
 -1.34082577e+10, -8.13188141e+09,  1.71687775e+10,
  4.22811127e+09,  4.34842606e+09,  2.43174597e+10,
  1.29962562e+11,  2.06104685e+10,  3.40849195e+08,
 -7.41704173e+10,  1.18621523e+10, -2.54342240e+10,
 -9.65473646e+09, -3.07997804e+10,  1.56097731e+10,
  1.05556369e+09, -1.06939824e+10, -3.88781845e+09,
 -8.14934114e+08,  4.31742313e+08,  2.58467057e+09,
  6.49520498e+09,  1.22026428e+10, -6.02102489e+09,
  4.57115666e+09, -9.06830161e+09,  2.44216691e+08,
 -9.19637628e+09, -1.89591728e+10,  4.09688141e+09,
 -8.13188141e+09, -4.31508638e+10,  3.70394070e+10,
  4.34842606e+09,  1.98506003e+10, -1.62421303e+10,
  2.06104685e+10,  7.66246120e+10, -1.07623426e+10,
  1.18621523e+10, -7.22475110e+09, -3.49576135e+09,
 -3.07997804e+10, -2.00125256e+10, -3.55609810e+09,
 -2.03524158e+10, -3.88781845e+09, -2.90717277e+10,
 -1.98956808e+09,  2.58467057e+09,  1.51200586e+09,
 -5.92587848e+09, -6.02102489e+09,  1.34434627e+10,
 -1.39335286e+10,  2.44216691e+08, -3.98836750e+10,
  1.01987555e+10,  4.09688141e+09, -1.93946779e+10,
  1.71687775e+10,  3.70394070e+10, -7.94883895e+10,
  2.43174597e+10, -1.62421303e+10,  1.91043045e+10,
  3.40849195e+08, -1.07623426e+10,  1.33216486e+11,
 -2.54342240e+10, -3.49576135e+09, -1.58661420e+10,
  1.56097731e+10, -3.55609810e+09,  1.64283529e+10,
  1.89209078e+10, -4.53257687e+09,  1.21753294e+10,
 -1.95445461e+10,  3.15871157e+09,  1.55626647e+10,
  8.46395315e+09, -3.12170571e+09,  2.14659590e+09,
 -1.11098463e+10, -2.76237189e+09,  2.03244510e+10,
 -2.11084114e+10, -1.16195095e+10,  1.88193265e+10,
 -3.29346413e+10,  9.27946559e+09, -1.70012811e+10,
  7.56508053e+08,  4.33459004e+09, -2.59982430e+10,
 -7.41704173e+10,  1.18621523e+10, -2.54342240e+10,
  1.41650000e+11, -3.29000000e+10,  1.77000000e+10,
 -1.09235066e+10,  2.63012445e+10, -1.82946193e+10,
 -4.53257687e+09,  5.28963397e+09, -2.47115666e+08,
  3.15871157e+09, -8.40183016e+08, -1.07735578e+10,
 -3.12170571e+09,  1.84171157e+10, -1.08464714e+10,
 -2.76237189e+09,  1.39227672e+09, -9.28814056e+09,
 -1.16195095e+10, -5.44731406e+10,  4.31375476e+10,
  9.27946559e+09, -2.95963543e+10,  1.40676720e+10,
  4.33459004e+09,  7.08971449e+09, -9.18295022e+09,
  1.18621523e+10, -7.22475110e+09, -3.49576135e+09,
 -3.29000000e+10,  1.17240000e+11, -3.94200000e+10,
  2.63012445e+10, -5.72943119e+10,  2.60487775e+10,
  1.21753294e+10, -2.47115666e+08,  1.73193265e+10,
  1.55626647e+10, -1.07735578e+10,  6.12466325e+09,
  2.14659590e+09, -1.08464714e+10,  2.22102928e+10,
  2.03244510e+10, -9.28814056e+09, -5.70721083e+09,
  1.88193265e+10,  4.31375476e+10, -9.78411786e+10,
 -1.70012811e+10,  1.40676720e+10, -4.85184919e+10,
 -2.59982430e+10, -9.18295022e+09, -3.29619253e+10,
 -2.54342240e+10, -3.49576135e+09, -1.58661420e+10,
  1.77000000e+10, -3.94200000e+10,  1.79360000e+11,
 -1.82946193e+10,  2.60487775e+10, -2.41193338e+10,
  1.28695461e+09,  5.66412884e+09,  1.20373353e+09,
  2.52797218e+09, -3.90585652e+09, -8.35614934e+08,
 -3.54082943e+10, -6.92774524e+09,  8.94033675e+09,
 -3.22868594e+10, -5.55988287e+09,  9.68521230e+09,
 -3.76033675e+09,  5.35805271e+09,  2.50583455e+10,
  7.14121523e+09,  6.42679356e+08, -2.20257687e+10,
 -4.72005124e+09, -1.67920937e+09, -3.49816984e+08,
 -9.65473646e+09, -3.07997804e+10,  1.56097731e+10,
 -1.09235066e+10,  2.63012445e+10, -1.82946193e+10,
  8.57976428e+10,  1.09063690e+10, -1.89915813e+10,
  5.66412884e+09,  4.75251830e+09,  2.99935578e+09,
 -3.90585652e+09,  6.62234993e+09, -4.50717423e+08,
 -6.92774524e+09, -5.43709297e+10,  3.46387262e+10,
 -5.55988287e+09, -4.06863470e+10,  2.77994143e+10,
  5.35805271e+09,  2.04473939e+10, -1.96452635e+10,
  6.42679356e+08,  7.28418741e+09, -1.03583968e+10,
 -1.67920937e+09, -1.39840922e+10, -2.94395315e+09,
 -3.07997804e+10, -2.00125256e+10, -3.55609810e+09,
  2.63012445e+10, -5.72943119e+10,  2.60487775e+10,
  1.09063690e+10,  1.47241757e+11, -5.45318448e+10,
  1.20373353e+09,  2.99935578e+09,  3.04903367e+09,
 -8.35614934e+08, -4.50717423e+08,  1.72892387e+09,
  8.94033675e+09,  3.46387262e+10, -7.83219107e+10,
  9.68521230e+09,  2.77994143e+10, -7.87758785e+10,
  2.50583455e+10, -1.96452635e+10,  1.66196047e+10,
 -2.20257687e+10, -1.03583968e+10, -2.77950952e+10,
 -3.49816984e+08, -2.94395315e+09, -5.77092972e+09,
  1.56097731e+10, -3.55609810e+09,  1.64283529e+10,
 -1.82946193e+10,  2.60487775e+10, -2.41193338e+10,
 -1.89915813e+10, -5.45318448e+10,  1.76957233e+11,
};

pylith::feassemble::ElasticityImplicitGravData3DQuadratic::ElasticityImplicitGravData3DQuadratic(void)
{ // constructor
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  matType = const_cast<char*>(_matType);
  matDBFilename = const_cast<char*>(_matDBFilename);
  matId = _matId;
  matLabel = const_cast<char*>(_matLabel);
  dt = _dt;
  gravityVec = const_cast<PylithScalar*>(_gravityVec);
  vertices = const_cast<PylithScalar*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDerivRef = const_cast<PylithScalar*>(_basisDerivRef);
  fieldTIncr = const_cast<PylithScalar*>(_fieldTIncr);
  fieldT = const_cast<PylithScalar*>(_fieldT);
  fieldTmdt = const_cast<PylithScalar*>(_fieldTmdt);
  valsResidual = const_cast<PylithScalar*>(_valsResidual);
  valsJacobian = const_cast<PylithScalar*>(_valsJacobian);
} // constructor

pylith::feassemble::ElasticityImplicitGravData3DQuadratic::~ElasticityImplicitGravData3DQuadratic(void)
{}


// End of file
