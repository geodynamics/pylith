// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#if !defined(pylith_bc_neumanndataquad4_hh)
#define pylith_bc_neumanndataquad4_hh

#include "NeumannData.hh"

namespace pylith {
  namespace bc {
     class NeumannDataQuad4;
  } // pylith
} // bc

class pylith::bc::NeumannDataQuad4 : public NeumannData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  NeumannDataQuad4(void);

  /// Destructor
  ~NeumannDataQuad4(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const char* _meshFilename;

  // Quadrature information
  static const int _numBasis;
  static const int _numQuadPts;
  static const PylithScalar _quadPts[];
  static const PylithScalar _quadWts[];
  static const PylithScalar _basis[];
  static const PylithScalar _basisDerivRef[];

  // BC information
  static const char* _spatialDBFilename;
  static const int _id;
  static const char* _label;

  // Mesh information
  static const int _spaceDim;
  static const int _cellDim;
  static const int _numVertices;
  static const int _numCells;
  static const int _numCorners;
  static const int _cells[];

  // Calculated values.
  static const PylithScalar _tractionsCell[];
  static const PylithScalar _valsResidual[];

};

#endif // pylith_bc_neumanndataquad4_hh

// End of file
