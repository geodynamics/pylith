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

#if !defined(pylith_bc_dirichletdatamultitri3_hh)
#define pylith_bc_dirichletdatamultitri3_hh

#include "DirichletDataMulti.hh"

namespace pylith {
  namespace bc {
     class DirichletDataMultiTri3;
  } // pylith
} // bc

class pylith::bc::DirichletDataMultiTri3 : public DirichletDataMulti
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  DirichletDataMultiTri3(void);

  /// Destructor
  ~DirichletDataMultiTri3(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const int _numDOF; ///< Number of degrees of freedom at each point.

  static const int _numFixedDOFA; ///< Number of fixedDOF at constrained points.
  static const int _numConstrainedPtsA; ///< Number of points constrained.
  static const int _idA; ///< Boundary condition identifier
  static const char* _labelA; /// Label for boundary condition group
  static const int _fixedDOFA[]; ///< Degrees of freedom constrained at points
  static const int _constrainedPointsA[]; ///< Array of indices of constrained pts.
  static const double _valuesA[]; ///< Values at constrained points.
  static const char* _dbFilenameA; ///< Filename of simple spatial database.

  static const int _numFixedDOFB; ///< Number of fixedDOF at constrained points.
  static const int _numConstrainedPtsB; ///< Number of points constrained.
  static const int _idB; ///< Boundary condition identifier
  static const char* _labelB; /// Label for boundary condition group
  static const int _fixedDOFB[]; ///< Degrees of freedom constrained at points
  static const int _constrainedPointsB[]; ///< Array of indices of constrained pts.
  static const double _valuesB[]; ///< Values at constrained points.
  static const char* _dbFilenameB; ///< Filename of simple spatial database.

  static const double _field[]; ///< Values in field
  static const int _constraintSizes[]; ///< Number of constrained DOF at each vertex
  static const int _constrainedDOF[]; ///< Indices of constrained DOF at each constrained vertex

  static const char* _meshFilename; ///< Filename of input mesh.
};

#endif // pylith_bc_dirichletdatamultitri3_hh

// End of file
