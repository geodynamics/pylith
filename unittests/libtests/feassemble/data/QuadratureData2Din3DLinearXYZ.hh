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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application quadratureapp.

#if !defined(pylith_feassemble_quadraturedata2din3dlinearxyz_hh)
#define pylith_feassemble_quadraturedata2din3dlinearxyz_hh

#include "QuadratureData.hh"

namespace pylith {
  namespace feassemble {
     class QuadratureData2Din3DLinearXYZ;
  } // pylith
} // feassemble

class pylith::feassemble::QuadratureData2Din3DLinearXYZ : public QuadratureData
{

public: 

  /// Constructor
  QuadratureData2Din3DLinearXYZ(void);

  /// Destructor
  ~QuadratureData2Din3DLinearXYZ(void);

private:

  static const int _numVertices;

  static const int _spaceDim;

  static const int _numCells;

  static const int _cellDim;

  static const int _numBasis;

  static const int _numQuadPts;

  static const PylithScalar _vertices[];

  static const int _cells[];

  static const PylithScalar _verticesRef[];

  static const PylithScalar _quadPtsRef[];

  static const PylithScalar _quadWts[];

  static const PylithScalar _quadPts[];

  static const PylithScalar _basis[];

  static const PylithScalar _basisDerivRef[];

  static const PylithScalar _basisDeriv[];

  static const PylithScalar _jacobian[];

  static const PylithScalar _jacobianDet[];

  static const PylithScalar _jacobianInv[];

};

#endif // pylith_feassemble_quadraturedata2din3dlinearxyz_hh

// End of file
