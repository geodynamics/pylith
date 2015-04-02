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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_bc_boundarymeshdatatri3_hh)
#define pylith_bc_boundarymeshdatatri3_hh

#include "BoundaryMeshData.hh"

namespace pylith {
  namespace bc {
     class BoundaryMeshDataTri3;
  } // pylith
} // bc

class pylith::bc::BoundaryMeshDataTri3 : public BoundaryMeshData
{

public: 

  /// Constructor
  BoundaryMeshDataTri3(void);

  /// Destructor
  ~BoundaryMeshDataTri3(void);

private:

  static const char* _filename;

  static const char* _bcLabel;

  static const char* _faultLabel;
  static const int _faultId;

  static const int _numCorners;
  static const int _numCells;
  static const bool _isSimplexMesh;

  static const int _numVerticesNoFault;

  static const int _numVerticesFault;

};

#endif // pylith_bc_boundarymeshdatatri3_hh

// End of file
