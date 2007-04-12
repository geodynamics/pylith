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
// This file was generated from python application integratorinertia3dlinear.

#if !defined(pylith_feassemble_integratordatainertia3dlinear_hh)
#define pylith_feassemble_integratordatainertia3dlinear_hh

#include "IntegratorData.hh"

namespace pylith {
  namespace feassemble {
     class IntegratorDataInertia3DLinear;
  } // pylith
} // feassemble

class pylith::feassemble::IntegratorDataInertia3DLinear : public IntegratorData
{

public: 

  /// Constructor
  IntegratorDataInertia3DLinear(void);

  /// Destructor
  ~IntegratorDataInertia3DLinear(void);

private:

  static const int _numVertices;

  static const int _spaceDim;

  static const int _numCells;

  static const int _cellDim;

  static const int _numBasis;

  static const int _numQuadPts;

  static const int _fiberDim;

  static const double _vertices[];

  static const int _cells[];

  static const double _quadPts[];

  static const double _quadWts[];

  static const double _basis[];

  static const double _basisDeriv[];

  static const double _fieldIn[];

  static const double _valsAction[];

  static const double _valsMatrix[];

  static const double _valsLumped[];

};

#endif // pylith_feassemble_integratordatainertia3dlinear_hh

// End of file
