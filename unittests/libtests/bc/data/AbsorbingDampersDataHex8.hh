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
  static const double _quadPts[];
  static const double _quadWts[];
  static const double _basis[];
  static const double _basisDerivRef[];

  static const char* _spatialDBFilename;
  static const int _id;
  static const char* _label;

  static const double _dt;
  static const double _fieldTIncr[];
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

#endif // pylith_bc_absorbingdampersdatahex8_hh

// End of file
