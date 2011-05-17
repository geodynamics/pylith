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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/*
  domain/cells (/topology/cells) [ncells, ncorners]
  domain/vertices (/geometry/vertices) [nvertices, spacedim]

  time (uniform time step)
    start, step, ntimesteps

  Repeat for each time slice
    topology (cell type, ncells, /Xdmf/Domain/DataItem[@Name="cells"])
    geometry (type, /Xdmf/Domain/DataItem[@Name="vertices"]
    Repeat for each field
      attribute (name, scalar, center, ntimesteps, npts, fiberdim)
      if 2-D vector, break up into components

*/

#if !defined(pylith_meshio_xdmf_hh)
#define pylith_meshio_xdmf_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include <fstream> // HASA std::ofstream
#include <string> // USES std::string

// Xdmf -----------------------------------------------------------------
/// Write Xdmf file associated with HDF5 file for use with VTK.
class pylith::meshio::Xdmf
{ // HDF5
  friend class TestXdmf; // Unit testing

// PUBLIC METHODS -------------------------------------------------------
public :

  /// Default constructor.
  Xdmf(void);

  /// Destructor
  ~Xdmf(void);

  /** Write Xdmf file associated with HDF5 file.
   *
   * @param filenameXdmf Name of Xdmf file.
   * @param filenameHDF5 Name of HDF5 file.
   */
  void write(const char* filenameXdfm,
	     const char* filenameHDF5);

// PRIVATE METHODS ------------------------------------------------------
private :

  /** Write domain cell information.
   *
   * @param numCells Number of cells.
   * @param numCorners Number of vertices in a cell.
   */
  void _writeDomainCells(const int numcells,
			 const int numCorners);

  /** Write domain vertices information.
   *
   * @param numVertices Number of vertices.
   * @param spaceDim Spatial dimension.
   */
  void _writeDomainVertices(const int numVertices,
			    const int spaceDim);

  /** Write grid topology information.
   *
   * @param cellType Name for cell type.
   * @param numCells Number of cells.
   */
  void _writeGridTopology(const char* cellType,
			  const int numells);

  /** Write Grid geometry.
   *
   * @param spaceDim Spatial dimension.
   */
  void _writeGridGeometry(const int spaceDim);

  /** Write grid attribute.
   * 
   * @param name Name of attribute.
   * @param center Vertex/Cell center.
   * @param numTimeSteps Number of time steps.
   * @param numPoints Number of vertices or cells.
   * @param fiberDim Fiber dimension for attribute.
   * @param iTime Index of time step.
   */
  void _writeGridAttribute(const char* name,
			   const char* center,
			   const int numTimeSteps,
			   const int numPoints,
			   const int fiberDim,
			   const int iTime);

// PRIVATE MEMBERS ------------------------------------------------------
private :

  std::ofstream _file; ///< Xdmf file

}; // Xdmf

#endif // pylith_meshio_xdmf_hh


// End of file 
