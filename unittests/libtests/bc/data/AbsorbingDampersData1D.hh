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

#if !defined(pylith_bc_absorbingdampersdata1d_hh)
#define pylith_bc_absorbingdampersdata1d_hh

#include "AbsorbingDampersData.hh"

namespace pylith {
  namespace bc {
     class AbsorbingDampersData1D;
  } // pylith
} // bc

class pylith::bc::AbsorbingDampersData1D : public AbsorbingDampersData
{

public: 

  /// Constructor
  AbsorbingDampersData1D(void);

  /// Destructor
  ~AbsorbingDampersData1D(void);

private:

  static const char* _meshFilename;

  static const int _spaceDim;

  static const int _cellDim;

  static const int _numVertices;

  static const int _numCells;

  static const int _numBasis;

  static const int _numQuadPts;

  static const char* _spatialDBFilename;

  static const double _dt;

  static const double _vertices[];

  static const int _cells[];

  static const double _verticesRef[];

  static const double _quadPts[];

  static const double _quadWts[];

  static const double _basis[];

  static const double _basisDerivRef[];

  static const double _fieldTpdt[];

  static const double _fieldT[];

  static const double _fieldTmdt[];

  static const double _valsResidual[];

  static const double _valsJacobian[];

};

#endif // pylith_bc_absorbingdampersdata1d_hh

// End of file
