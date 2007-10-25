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

#if !defined(pylith_bc_absorbingdampersdataquad4_hh)
#define pylith_bc_absorbingdampersdataquad4_hh

#include "AbsorbingDampersData.hh"

namespace pylith {
  namespace bc {
     class AbsorbingDampersDataQuad4;
  } // pylith
} // bc

class pylith::bc::AbsorbingDampersDataQuad4 : public AbsorbingDampersData
{

public: 

  /// Constructor
  AbsorbingDampersDataQuad4(void);

  /// Destructor
  ~AbsorbingDampersDataQuad4(void);

private:

  static const char* _meshFilename;

  static const int _numBasis;
  static const int _numQuadPts;
  static const double _quadPts[];
  static const double _quadWts[];
  static const double _basis[];
  static const double _basisDerivRef[];

  static const char* _spatialDBFilename;
  static const int _id;
  static const char* _label;

  static const double _dt;
  static const double _fieldTpdt[];
  static const double _fieldT[];
  static const double _fieldTmdt[];

  static const int _spaceDim;
  static const int _cellDim;
  static const int _numVertices;
  static const int _numCells;
  static const int _numCorners;
  static const int _cells[];

  static const double _dampingConsts[];
  static const double _valsResidual[];
  static const double _valsJacobian[];

};

#endif // pylith_bc_absorbingdampersdataquad4_hh

// End of file
