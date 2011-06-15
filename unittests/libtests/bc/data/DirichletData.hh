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

#if !defined(pylith_bc_dirichletdata_hh)
#define pylith_bc_dirichletdata_hh

namespace pylith {
  namespace bc {
     class DirichletData;
  } // pylith
} // bc

class pylith::bc::DirichletData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  DirichletData(void);

  /// Destructor
  ~DirichletData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  double tRef; ///< Reference time for rate of change of values
  double valueRate; ///< Rate of change of value at constrained points.

  int numDOF; ///< Number of degrees of freedom at each point.
  int numFixedDOF; ///< Number of fixedDOF at constrained points.
  int numConstrainedPts; ///< Number of points constrained.

  int id; ///< Boundary condition identifier
  char* label; ///< Label for boundary condition group

  int* fixedDOF; ///< Degrees of freedom that are constrained at each point
  int* constrainedPoints; ///< Array of indices of constrained points.
  double* valuesInitial; ///< Values at constrained points.

  char* meshFilename; ///< Filename for input mesh.
  char* dbFilename; ///< Filename of simple spatial database.
};

#endif // pylith_bc_dirichletdata_hh

// End of file
