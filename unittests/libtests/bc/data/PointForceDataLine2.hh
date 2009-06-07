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

#if !defined(pylith_bc_forcedatadataline2_hh)
#define pylith_bc_forcedatadataline2_hh

#include "PointForceData.hh"

namespace pylith {
  namespace bc {
     class PointForceDataLine2;
  } // pylith
} // bc

class pylith::bc::PointForceDataLine2 : public PointForceData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  PointForceDataLine2(void);

  /// Destructor
  ~PointForceDataLine2(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const double _tRef; ///< Reference time for rate of change of forces.
  static const double _forceRate; ///< Rate of change of force.
  static const double _tResidual; ///< Time for computing residual.

  static const int _numDOF; ///< Number of degrees of freedom at each point.
  static const int _numForceDOF; ///< Number of forces at points.
  static const int _numForcePts; ///< Number of points with forces.

  static const int _id; ///< Boundary condition identifier
  static const char* _label; ///< Label for boundary condition group

  static const int _forceDOF[]; ///< Degrees of freedom that are constrained at each point
  static const int _forcePoints[]; ///< Array of indices of points with forces.
  static const double _forceInitial[]; ///< Forces at points.
  static const double _residual[]; ///< Residual field.

  static const char* _meshFilename; ///< Filename for input mesh.
  static const char* _dbFilename; ///< Filename of simple spatial database.
};

#endif // pylith_bc_forcedatadataline2_hh

// End of file
