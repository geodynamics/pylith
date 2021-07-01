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

#if !defined(pylith_bc_absorbingdampersdatahex8_hh)
#define pylith_bc_absorbingdampersdatahex8_hh

#include "AbsorbingDampersData.hh"

namespace pylith {
  namespace bc {
     class AbsorbingDampersDataHex8;
  } // pylith
} // bc

class pylith::bc::AbsorbingDampersDataHex8 : public AbsorbingDampersData
{

public: 

  /// Constructor
  AbsorbingDampersDataHex8(void);

  /// Destructor
  ~AbsorbingDampersDataHex8(void);

private:

  static const char* _meshFilename;

  static const int _numBasis;
  static const int _numQuadPts;
  static const PylithScalar _quadPts[];
  static const PylithScalar _quadWts[];
  static const PylithScalar _basis[];
  static const PylithScalar _basisDerivRef[];

  static const char* _spatialDBFilename;
  static const int _id;
  static const char* _label;

  static const PylithScalar _dt;
  static const PylithScalar _fieldTIncr[];
  static const PylithScalar _fieldT[];
  static const PylithScalar _fieldTmdt[];

  static const int _spaceDim;
  static const int _cellDim;
  static const int _numVertices;
  static const int _numCells;
  static const int _numCorners;
  static const int _cells[];

  static const PylithScalar _dampingConsts[];
  static const PylithScalar _valsResidual[];
  static const PylithScalar _valsJacobian[];

};

#endif // pylith_bc_absorbingdampersdatahex8_hh

// End of file
