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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_faults_cohesiveimpulsesdata_hh)
#define pylith_faults_cohesiveimpulsesdata_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
  namespace faults {
     class CohesiveImpulsesData;
  } // pylith
} // faults

class pylith::faults::CohesiveImpulsesData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  CohesiveImpulsesData(void);

  /// Destructor
  ~CohesiveImpulsesData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* meshFilename; ///< Filename for input mesh

  /// @name Quadrature information
  //@{
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int cellDim; ///< Number of dimensions associated with cell
  int numBasis; ///< Number of vertices in cell
  int numQuadPts; ///< Number of quadrature points
  PylithScalar* quadPts; ///< Coordinates of quad pts in ref cell
  PylithScalar* quadWts; ///< Weights of quadrature points
  PylithScalar* basis; ///< Basis fns at quadrature points
  PylithScalar* basisDeriv; ///< Derivatives of basis fns at quad pts
  PylithScalar* verticesRef; ///< Coordinates of vertices in ref cell (dual basis)
  //@}

  /// @name Fault information
  //@{
  int id; ///< Fault material identifier
  char* label; ///< Label for fault
  char* impulseAmpFilename; ///< Name of db for final slip
  int* impulseDOF; ///< Fault DOF associated with impulses.
  int numComponents; ///< Number of components in impulses.
  //@}

  /// @name Input fields
  //@{
  PylithScalar* fieldT; ///< Field over domain at time t.
  PylithScalar* fieldIncr; ///< Solution increment field over domain at time t.
  //@}

  /// @name Calculated values.
  //@{
  PylithScalar* orientation; ///< Expected values for fault orientation.
  PylithScalar* area; ///< Expected values for fault area.
  PylithScalar* amplitude; ///< Expected values for impulse amplitude.
  int numImpulses; ///< Number of impulses.

  /// Expected values from residual calculation using solution increment.
  PylithScalar* residualIncr;

  int* constraintVertices; ///< Expected points for constraint vertices
  int numConstraintVert; ///< Number of constraint vertices
  //@}

};

#endif // pylith_faults_cohesiveimpulsesdata_hh

// End of file
