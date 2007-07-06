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

#if !defined(pylith_faults_boundarydatahex8_hh)
#define pylith_faults_boundarydatahex8_hh

#include "BoundaryData.hh"

namespace pylith {
  namespace faults {
     class BoundaryDataHex8;
  } // pylith
} // faults

class pylith::faults::BoundaryDataHex8 : public BoundaryData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  BoundaryDataHex8(void);

  /// Destructor
  ~BoundaryDataHex8(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const int _numVertices; ///< Number of vertices
  static const int _spaceDim; ///< Number of dimensions in vertex coordinates
  static const int _numCells; ///< Number of cells
  static const int _cellDim; ///< Number of dimensions associated with cell

  static const double _vertices[]; ///< Pointer to coordinates of vertices
  static const int _numCorners[]; ///< Number of vertices in cell
  static const int _cells[]; ///< Pointer to indices of vertices in cells

  static const char* _filename; ///< Filename of input mesh
};

#endif // pylith_faults_boundarydatahex8_hh

// End of file
