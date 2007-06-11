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

#if !defined(pylith_faults_cohesivedatatri3d_hh)
#define pylith_faults_cohesivedatatri3d_hh

#include "CohesiveData.hh"

namespace pylith {
  namespace faults {
     class CohesiveDataTri3d;
  } // pylith
} // faults

class pylith::faults::CohesiveDataTri3d : public CohesiveData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  CohesiveDataTri3d(void);

  /// Destructor
  ~CohesiveDataTri3d(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const int _numVertices; ///< Number of vertices
  static const int _spaceDim; ///< Number of dimensions in vertex coordinates
  static const int _numCells; ///< Number of cells
  static const int _cellDim; ///< Number of dimensions associated with cell

  static const double _vertices[]; ///< Pointer to coordinates of vertices
  static const int _numCorners[]; ///< Number of vertices in cell
  static const int _cells[]; ///< Pointer to indices of vertices in cells
  static const int _materialIds[]; ///< Pointer to cell material identifiers

  static const int _groups[]; ///< Groups of points
  static const int _groupSizes[]; ///< Sizes of groups
  static const char* _groupNames[]; ///< Array of group names
  static const char* _groupTypes[]; ///< Array of group types
  static const int _numGroups; ///< Number of groups

  static const char* _filename; ///< Filename of input mesh
};

#endif // pylith_faults_cohesivedatatri3d_hh

// End of file
