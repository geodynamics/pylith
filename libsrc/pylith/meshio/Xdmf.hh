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
// Copyright (c) 2010-2015 University of California, Davis
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

#include "pylith/utils/arrayfwd.hh" // USES scalar_array, std::vector

#include <vector> // USES std::vector
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
  void write(const char* filenameXdmf,
	     const char* filenameHDF5);

// PRIVATE STRUCTS ------------------------------------------------------
public :

  /// Metadata associated with fields.
  struct FieldMetadata {
    std::string name; ///< Name of field.
    std::string vectorFieldType; ///< Type of field.
    std::string domain; ///< Domain over which field is given.
    int numPoints; ///< Number of points in field.
    int fiberDim; ///< Number of values in field at each point, time.
    int numTimeSteps; ///< Number of time steps for field.
  }; // FieldMetadata

// PRIVATE METHODS ------------------------------------------------------
private :

  /** Get timestamps from HDF5 file.
   *
   * @param timeStamps Array of time stamps.
   * @param h5 HDF5 file.
   */
  void _getTimeStamps(scalar_array* timeStamps,
		      HDF5& h5);

  /** Get field metadata from HDF5 file.
   *
   * @param metadata Array of field metadata.
   * @param h5 HDF5 file.
   */
  void _getFieldMetadata(std::vector<FieldMetadata>* metadata,
			 HDF5& h5);

  /** Write domain cell information.
   *
   * @param numCells Number of cells.
   * @param numCorners Number of vertices in a cell.
   */
  void _writeDomainCells(const int numCells,
			 const int numCorners);

  /** Write domain vertices information.
   *
   * @param numVertices Number of vertices.
   * @param spaceDim Spatial dimension.
   */
  void _writeDomainVertices(const int numVertices,
			    const int spaceDim);

  /** Write time stamps associated with fields.
   *
   * @param timeStamps Array of time stamps.
   */
  void _writeTimeStamps(const scalar_array& timeStamps);

  /** Write grid topology information.
   *
   * @param cellType Name for cell type.
   * @param numCells Number of cells.
   */
  void _writeGridTopology(const char* cellType,
			  const int numCells);

  /** Write Grid geometry.
   *
   * @param spaceDim Spatial dimension.
   */
  void _writeGridGeometry(const int spaceDim);

  /** Write grid attribute.
   * 
   * @param metadata Metadata for field.
   * @param iTime Index of time step.
   */
  void _writeGridAttribute(const FieldMetadata& metadata,
			   const int iTime);

  /** Write grid attribute as single component (for 2-D vector).
   * 
   * @param metadata Metadata for field.
   * @param iTime Index of time step.
   * @param component Index of component.
   * @param spaceDim Spatial dimension.
   */
  void _writeGridAttributeComponent(const FieldMetadata& metadata,
				    const int iTime,
				    const int component,
				    const int spaceDim);

// PRIVATE MEMBERS ------------------------------------------------------
private :

  std::ofstream _file; ///< Xdmf file

}; // Xdmf

#endif // pylith_meshio_xdmf_hh


// End of file 
