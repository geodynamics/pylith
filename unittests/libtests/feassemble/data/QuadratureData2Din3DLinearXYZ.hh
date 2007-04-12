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
// This file was generated from python application quadrature2din3dlinear.

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

  static const double _vertices[];

  static const int _cells[];

  static const double _quadPtsRef[];

  static const double _quadWts[];

  static const double _basis[];

  static const double _basisDeriv[];

  static const double _quadPts[];

  static const double _jacobian[];

  static const double _jacobianDet[];

  static const double _jacobianInv[];

};

#endif // pylith_feassemble_quadraturedata2din3dlinearxyz_hh

// End of file
