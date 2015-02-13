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

#if !defined(pylith_bc_boundarymeshdataquad4_hh)
#define pylith_bc_boundarymeshdataquad4_hh

#include "BoundaryMeshData.hh"

namespace pylith {
  namespace bc {
     class BoundaryMeshDataQuad4;
  } // pylith
} // bc

class pylith::bc::BoundaryMeshDataQuad4 : public BoundaryMeshData
{

public: 

  /// Constructor
  BoundaryMeshDataQuad4(void);

  /// Destructor
  ~BoundaryMeshDataQuad4(void);

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

#endif // pylith_bc_boundarymeshdataquad4_hh

// End of file
