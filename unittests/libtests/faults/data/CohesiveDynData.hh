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

#if !defined(pylith_faults_cohesivedynldata_hh)
#define pylith_faults_cohesivedynldata_hh

namespace pylith {
  namespace faults {
     class CohesiveDynData;
  } // pylith
} // faults

class pylith::faults::CohesiveDynData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  CohesiveDynData(void);

  /// Destructor
  ~CohesiveDynData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* meshFilename; ///< Filename for input mesh

  /// @name Quadrature information
  //@{
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int cellDim; ///< Number of dimensions associated with cell
  int numBasis; ///< Number of vertices in cell
  int numQuadPts; ///< Number of quadrature points
  double* quadPts; ///< Coordinates of quad pts in ref cell
  double* quadWts; ///< Weights of quadrature points
  double* basis; ///< Basis fns at quadrature points
  double* basisDeriv; ///< Derivatives of basis fns at quad pts
  double* verticesRef; ///< Coordinates of vertices in ref cell (dual basis)
  //@}

  /// @name Fault information
  //@{
  int id; ///< Fault material identifier
  char* label; ///< Label for fault
  char* initialTractFilename; ///< Name of db for initial tractions.
  //@}

  /// @name Input fields
  //@{
  double* fieldT; ///< Solution field at time t.
  double* fieldIncrStick; ///< Soln increment field at time t for stick case.
  double* fieldIncrSlip; ///< Soln increment field at time t for slipping case.
  double* fieldIncrOpen; ///< Soln increment field at time t for opening case.
  double* jacobian; ///< Jacobian sparse matrix.
  //@}

  /// @name Calculated values.
  //@{
  double* orientation; ///< Expected values for fault orientation.
  double* area; ///< Expected values for fault area.
  double* forcesInitial; ///< Expected values for initial forces.
  double* fieldIncrSlipE; ///< Expected values for solution increment for slipping case.
  double* slipSlipE; ///< Expected values for slip for slipping case.
  double* fieldIncrOpenE; ///< Expected values for solution increment for opening case.
  double* slipOpenE; ///< Expected values for slip for opening case.

  int* constraintVertices; ///< Expected points for constraint vertices
  int numConstraintVert; ///< Number of constraint vertices
  //@}

};

#endif // pylith_faults_cohesivedynldata_hh

// End of file
