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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_faults_cohesivedata_hh)
#define pylith_faults_cohesivedata_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
  namespace faults {
     class CohesiveData;
  } // pylith
} // faults

class pylith::faults::CohesiveData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  CohesiveData(void);

  /// Destructor
  ~CohesiveData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  int numVertices; ///< Number of vertices
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int numCells; ///< Number of cells
  int cellDim; ///< Number of dimensions associated with cell

  PylithScalar* vertices; ///< Pointer to coordinates of vertices
  int* numCorners; ///< Number of vertices in cell
  int* cells; ///< Pointer to indices of vertices in cells
  int* materialIds; ///< Pointer to cell material identifiers

  int* groups; ///< Array of pointers to indices of points in groups
  int* groupSizes; ///< Array of sizes of each group
  char** groupNames; ///< Array of group names
  char** groupTypes; ///< Array of group types
  int numGroups; ///< Number of groups

  char* filename; ///< Filename for input mesh
};

#endif // pylith_faults_cohesivedata_hh

// End of file
