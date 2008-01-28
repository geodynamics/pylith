// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#if !defined(pylith_bc_dirichletpointsdataline2_hh)
#define pylith_bc_dirichletpointsdataline2_hh

#include "DirichletPointsData.hh"

namespace pylith {
  namespace bc {
     class DirichletPointsDataLine2;
  } // pylith
} // bc

class pylith::bc::DirichletPointsDataLine2 : public DirichletPointsData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  DirichletPointsDataLine2(void);

  /// Destructor
  ~DirichletPointsDataLine2(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const int _numDOF; ///< Number of degrees of freedom at each point.

  static const int _numFixedDOF; ///< Number of fixedDOF at constrained points.
  static const int _numConstrainedPts; ///< Number of points constrained.

  static const int _id; ///< Boundary condition identifier
  static const char* _label; /// Label for boundary condition group

  static const int _fixedDOF[]; ///< Degrees of freedom constrained at points

  static const int _constrainedPoints[]; ///< Array of indices of constrained pts.
  static const double _values[]; ///< Values at constrained points.

  static const char* _meshFilename; ///< Filename of input mesh.
  static const char* _dbFilename; ///< Filename of simple spatial database.
};

#endif // pylith_bc_dirichletpointsdataline2_hh

// End of file
