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
// This file was generated from python application quadratureapp.

#if !defined(pylith_feassemble_quadraturedata1din2dlinear_hh)
#define pylith_feassemble_quadraturedata1din2dlinear_hh

#include "QuadratureData.hh"

namespace pylith {
  namespace feassemble {
     class QuadratureData1Din2DLinear;
  } // pylith
} // feassemble

class pylith::feassemble::QuadratureData1Din2DLinear : public QuadratureData
{

public: 

  /// Constructor
  QuadratureData1Din2DLinear(void);

  /// Destructor
  ~QuadratureData1Din2DLinear(void);

private:

  static const int _numVertices;

  static const int _spaceDim;

  static const int _numCells;

  static const int _cellDim;

  static const int _numBasis;

  static const int _numQuadPts;

  static const double _vertices[];

  static const int _cells[];

  static const double _verticesRef[];

  static const double _quadPtsRef[];

  static const double _quadWts[];

  static const double _quadPts[];

  static const double _basis[];

  static const double _basisDerivRef[];

  static const double _basisDeriv[];

  static const double _jacobian[];

  static const double _jacobianDet[];

  static const double _jacobianInv[];

};

#endif // pylith_feassemble_quadraturedata1din2dlinear_hh

// End of file
