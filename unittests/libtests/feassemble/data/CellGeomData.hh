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

#if !defined(pylith_feassemble_cellgeomdata_hh)
#define pylith_feassemble_cellgeomdata_hh

namespace pylith {
  namespace feassemble {
     class CellGeomData;
  } // pylith
} // feassemble

class pylith::feassemble::CellGeomData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  CellGeomData(void);

  /// Destructor
  ~CellGeomData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  int cellDim; ///< Number of dimensions associated with cell
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int numCorners; ///< Number of vertices in cell

  int numLocs; ///< Number of locations

  double* gravityVec; ///< Gravity vector for problem
  double* vertices; ///< Coordinates of vertices of cell
  double* locations; ///< Locations where Jacobian is computed
  double* jacobian; ///< Jacobian at locations
  double* jacobianDet; ///< Determinant of Jacobian at locations

};

#endif // pylith_feassemble_cellgeomdata_hh


// End of file
