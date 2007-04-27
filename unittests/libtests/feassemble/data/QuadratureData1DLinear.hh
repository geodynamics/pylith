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
// This file was generated from python application quadrature1dlinear.

#if !defined(pylith_feassemble_quadraturedata1dlinear_hh)
#define pylith_feassemble_quadraturedata1dlinear_hh

#include "QuadratureData.hh"

namespace pylith {
  namespace feassemble {
     class QuadratureData1DLinear;
  } // pylith
} // feassemble

class pylith::feassemble::QuadratureData1DLinear : public QuadratureData
{

public: 

  /// Constructor
  QuadratureData1DLinear(void);

  /// Destructor
  ~QuadratureData1DLinear(void);

private:

  static const int _numVertices;

  static const int _spaceDim;

  static const int _numCells;

  static const int _cellDim;

  static const int _numBasis;

  static const int _numQuadPts;

  static const double _vertices[];

  static const int _cells[];

  static const double _quadPtsRef[];

  static const double _quadWts[];

  static const double _quadPts[];

  static const double _basisVert[];

  static const double _basisDerivVert[];

  static const double _jacobianVert[];

  static const double _jacobianDetVert[];

  static const double _basisQuad[];

  static const double _basisDerivQuad[];

  static const double _jacobianQuad[];

  static const double _jacobianDetQuad[];

  static const double _jacobianInvQuad[];

};

#endif // pylith_feassemble_quadraturedata1dlinear_hh

// End of file
