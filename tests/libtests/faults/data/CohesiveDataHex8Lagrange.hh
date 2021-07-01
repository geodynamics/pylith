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

#if !defined(pylith_faults_cohesivedatahex8lagrange_hh)
#define pylith_faults_cohesivedatahex8lagrange_hh

#include "CohesiveData.hh"

namespace pylith {
  namespace faults {
     class CohesiveDataHex8Lagrange;
  } // pylith
} // faults

class pylith::faults::CohesiveDataHex8Lagrange : public CohesiveData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  CohesiveDataHex8Lagrange(void);

  /// Destructor
  ~CohesiveDataHex8Lagrange(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const int _numVertices; ///< Number of vertices
  static const int _spaceDim; ///< Number of dimensions in vertex coordinates
  static const int _numCells; ///< Number of cells
  static const int _cellDim; ///< Number of dimensions associated with cell

  static const int _numCorners[]; ///< Number of vertices in cell
  static const int _materialIds[]; ///< Pointer to cell material identifiers

  static const int _groupSizes[]; ///< Sizes of groups
  static const char* _groupNames[]; ///< Array of group names
  static const char* _groupTypes[]; ///< Array of group types
  static const int _numGroups; ///< Number of groups

  static const char* _filename; ///< Filename of input mesh
};

#endif // pylith_faults_cohesivedatahex8lagrange_hh

// End of file
