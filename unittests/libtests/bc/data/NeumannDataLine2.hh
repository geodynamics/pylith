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

#if !defined(pylith_bc_neumanndataline2_hh)
#define pylith_bc_neumanndataline2_hh

#include "NeumannData.hh"

namespace pylith {
  namespace bc {
     class NeumannDataLine2;
  } // pylith
} // bc

class pylith::bc::NeumannDataLine2 : public NeumannData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  NeumannDataLine2(void);

  /// Destructor
  ~NeumannDataLine2(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const char* _meshFilename;

  // Quadrature information
  static const int _numBasis;
  static const int _numQuadPts;
  static const double _quadPts[];
  static const double _quadWts[];
  static const double _basis[];
  static const double _basisDerivRef[];

  // BC information
  static const char* _spatialDBFilename;
  static const int _id;
  static const char* _label;

  // Mesh information
  static const int _spaceDim;
  static const int _cellDim;
  static const int _numBoundaryVertices;
  static const int _numBoundaryCells;
  static const int _numCorners;
  static const double _cellVertices[];

  // Calculated values.
  static const double _tractionsCell[];
  static const double _valsResidual[];

};

#endif // pylith_bc_neumanndataline2_hh

// End of file
