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

// DO NOT EDIT THIS FILE
// This file was generated from python application elasticityimplicitgrav.

#include "ElasticityImplicitGravData3DLinear.hh"

const int pylith::feassemble::ElasticityImplicitGravData3DLinear::_spaceDim = 3;

const int pylith::feassemble::ElasticityImplicitGravData3DLinear::_cellDim = 3;

const int pylith::feassemble::ElasticityImplicitGravData3DLinear::_numVertices = 4;

const int pylith::feassemble::ElasticityImplicitGravData3DLinear::_numCells = 1;

const int pylith::feassemble::ElasticityImplicitGravData3DLinear::_numBasis = 4;

const int pylith::feassemble::ElasticityImplicitGravData3DLinear::_numQuadPts = 1;

const char* pylith::feassemble::ElasticityImplicitGravData3DLinear::_matType = "ElasticIsotropic3D";

const char* pylith::feassemble::ElasticityImplicitGravData3DLinear::_matDBFilename = "data/elasticisotropic3d.spatialdb";

const int pylith::feassemble::ElasticityImplicitGravData3DLinear::_matId = 0;

const char* pylith::feassemble::ElasticityImplicitGravData3DLinear::_matLabel = "elastic isotropic 3-D";

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_dt =   1.00000000e-02;

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_gravityVec[] = {
  0.00000000e+00,  0.00000000e+00, -1.00000000e+08,
};

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_vertices[] = {
 -5.00000000e-01, -1.00000000e+00, -5.00000000e-01,
  2.00000000e+00, -5.00000000e-01, -4.00000000e-01,
  1.00000000e+00, -1.00000000e-01, -3.00000000e-01,
 -2.00000000e-01,  5.00000000e-01,  2.00000000e+00,
};

const int pylith::feassemble::ElasticityImplicitGravData3DLinear::_cells[] = {
0,1,2,3,
};

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_quadPts[] = {
 -5.00000000e-01, -5.00000000e-01, -5.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_quadWts[] = {
  1.33333333e+00,
};

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_basis[] = {
 -2.50000000e-01,  2.50000000e-01,  2.50000000e-01,
  2.50000000e-01,};

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_basisDerivRef[] = {
 -5.00000000e-01, -5.00000000e-01, -5.00000000e-01,
  5.00000000e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  5.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  5.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_fieldTIncr[] = {
  3.00000000e-01,  2.00000000e-01, -5.00000000e-01,
 -3.00000000e-01, -4.00000000e-01, -6.00000000e-01,
  2.00000000e-01,  6.00000000e-01,  3.00000000e-01,
 -6.00000000e-01, -1.00000000e-01, -3.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_fieldT[] = {
  8.00000000e-01,  1.00000000e-01, -6.00000000e-01,
 -1.00000000e-01, -2.00000000e-01, -5.00000000e-01,
  1.00000000e-01,  7.00000000e-01,  2.00000000e-01,
 -5.00000000e-01, -0.00000000e+00, -2.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_fieldTmdt[] = {
  1.00000000e-01,  1.00000000e-01, -3.00000000e-01,
 -2.00000000e-01, -1.00000000e-01, -5.00000000e-01,
  2.00000000e-01,  4.00000000e-01,  1.00000000e-01,
 -4.00000000e-01, -1.00000000e-01, -1.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_valsResidual[] = {
 -6.04851301e+09,  7.19421933e+10,  4.07639095e+10,
  1.11995353e+10,  1.19496190e+11,  2.47670074e+10,
 -1.62946097e+10, -1.94715799e+11, -1.17458953e+11,
  1.11435874e+10,  3.27741636e+09, -1.53219641e+10,
};

const double pylith::feassemble::ElasticityImplicitGravData3DLinear::_valsJacobian[] = {
  1.08203222e+10,  5.80793061e+09, -3.19702602e+08,
  5.03531599e+09, -1.66914498e+09, -2.27509294e+09,
 -1.72763321e+10, -4.57125155e+09,  4.59107807e+09,
  1.42069393e+09,  4.32465923e+08, -1.99628253e+09,
  5.80793061e+09,  2.32515489e+10, -8.10408922e+08,
 -1.66914498e+09,  3.51505576e+10, -7.37174721e+09,
 -4.57125155e+09, -6.22131351e+10,  1.33122677e+10,
  4.32465923e+08,  3.81102850e+09, -5.13011152e+09,
 -3.19702602e+08, -8.10408922e+08,  8.57372986e+09,
 -2.27509294e+09, -7.37174721e+09,  1.09665428e+10,
  4.59107807e+09,  1.33122677e+10, -2.15452292e+10,
 -1.99628253e+09, -5.13011152e+09,  2.00495663e+09,
  5.03531599e+09, -1.66914498e+09, -2.27509294e+09,
  4.48327138e+10, -2.22908922e+10,  1.19609665e+10,
 -5.65594796e+10,  2.50743494e+10, -1.42472119e+10,
  6.69144981e+09, -1.11431227e+09,  4.56133829e+09,
 -1.66914498e+09,  3.51505576e+10, -7.37174721e+09,
 -2.22908922e+10,  7.52342007e+10, -2.26338290e+10,
  2.50743494e+10, -1.21016729e+11,  3.96524164e+10,
 -1.11431227e+09,  1.06319703e+10, -9.64684015e+09,
 -2.27509294e+09, -7.37174721e+09,  1.09665428e+10,
  1.19609665e+10, -2.26338290e+10,  4.51979554e+10,
 -1.42472119e+10,  3.96524164e+10, -7.19962825e+10,
  4.56133829e+09, -9.64684015e+09,  1.58317844e+10,
 -1.72763321e+10, -4.57125155e+09,  4.59107807e+09,
 -5.65594796e+10,  2.50743494e+10, -1.42472119e+10,
  8.56232962e+10, -2.11957869e+10,  1.22676580e+10,
 -1.17874845e+10,  6.92688971e+08, -2.61152416e+09,
 -4.57125155e+09, -6.22131351e+10,  1.33122677e+10,
  2.50743494e+10, -1.21016729e+11,  3.96524164e+10,
 -2.11957869e+10,  2.01727385e+11, -6.93680297e+10,
  6.92688971e+08, -1.84975217e+10,  1.64033457e+10,
  4.59107807e+09,  1.33122677e+10, -2.15452292e+10,
 -1.42472119e+10,  3.96524164e+10, -7.19962825e+10,
  1.22676580e+10, -6.93680297e+10,  1.22023544e+11,
 -2.61152416e+09,  1.64033457e+10, -2.84820322e+10,
  1.42069393e+09,  4.32465923e+08, -1.99628253e+09,
  6.69144981e+09, -1.11431227e+09,  4.56133829e+09,
 -1.17874845e+10,  6.92688971e+08, -2.61152416e+09,
  3.67534077e+09, -1.08426270e+07,  4.64684015e+07,
  4.32465923e+08,  3.81102850e+09, -5.13011152e+09,
 -1.11431227e+09,  1.06319703e+10, -9.64684015e+09,
  6.92688971e+08, -1.84975217e+10,  1.64033457e+10,
 -1.08426270e+07,  4.05452292e+09, -1.62639405e+09,
 -1.99628253e+09, -5.13011152e+09,  2.00495663e+09,
  4.56133829e+09, -9.64684015e+09,  1.58317844e+10,
 -2.61152416e+09,  1.64033457e+10, -2.84820322e+10,
  4.64684015e+07, -1.62639405e+09,  1.06452912e+10,
};

pylith::feassemble::ElasticityImplicitGravData3DLinear::ElasticityImplicitGravData3DLinear(void)
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
  gravityVec = const_cast<double*>(_gravityVec);
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<double*>(_verticesRef);
  quadPts = const_cast<double*>(_quadPts);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDerivRef = const_cast<double*>(_basisDerivRef);
  fieldTIncr = const_cast<double*>(_fieldTIncr);
  fieldT = const_cast<double*>(_fieldT);
  fieldTmdt = const_cast<double*>(_fieldTmdt);
  valsResidual = const_cast<double*>(_valsResidual);
  valsJacobian = const_cast<double*>(_valsJacobian);
} // constructor

pylith::feassemble::ElasticityImplicitGravData3DLinear::~ElasticityImplicitGravData3DLinear(void)
{}


// End of file
