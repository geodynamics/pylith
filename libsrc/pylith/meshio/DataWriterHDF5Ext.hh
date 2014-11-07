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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/DataWriterHDF5Ext.hh
 *
 * @brief Object for writing finite-element data to HDF5 file with
 * datasets stored in raw external data files.
 *
 * HDF5 schema for PyLith output.
 *
 * / - root group
 *   geometry - group
       coordsys - attribute string with coordinate system
 *     vertices - dataset [nvertices, spacedim]
 *   topology - group
 *     cell_type - attribute string with cell type
 *     cells - dataset [ncells, ncorners]
 *   vertex_fields - group
 *     VERTEX_FIELD (name of vertex field) - dataset 
 *       [ntimesteps, nvertices, fiberdim]
 *   cell_fields - group
 *     CELL_FIELD (name of cell field) - dataset
 *       [ntimesteps, ncells, fiberdim]
 */

#if !defined(pylith_meshio_datawriterhdf5ext_hh)
#define pylith_meshio_datawriterhdf5ext_hh

// Include directives ---------------------------------------------------
#include "DataWriter.hh" // ISA DataWriter

#include <string> // USES std::string
#include <map> // HASA std::map

// DataWriterHDF5Ext ----------------------------------------------------
/// Object for writing finite-element data to HDF5 file.
class pylith::meshio::DataWriterHDF5Ext : public DataWriter
{ // DataWriterHDF5Ext
  friend class TestDataWriterHDF5ExtMesh; // unit testing
  friend class TestDataWriterHDF5ExtSubMesh; // unit testing
  friend class TestDataWriterHDF5ExtPoints; // unit testing
  friend class TestDataWriterHDF5ExtBCMesh; // unit testing
  friend class TestDataWriterHDF5ExtFaultMesh; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  DataWriterHDF5Ext(void);

  /// Destructor
  ~DataWriterHDF5Ext(void);

  /** Make copy of this object.
   *
   * @returns Copy of this.
   */
  DataWriter* clone(void) const;

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Set filename for HDF5 file.
   *
   * @param filename Name of HDF5 file.
   */
  void filename(const char* filename);

  /** Prepare for writing files.
   *
   * @param mesh Finite-element mesh. 
   * @param numTimeSteps Expected number of time steps for fields.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  void open(const topology::Mesh& mesh,
	    const int numTimeSteps,
	    const char* label =0,
	    const int labelId =0);

  /// Close output files.
  void close(void);

  /** Write field over vertices to file.
   *
   * @param t Time associated with field.
   * @param field Field over vertices.
   * @param mesh Mesh associated with output.
   */
  void writeVertexField(const PylithScalar t,
			topology::Field& field,
			const topology::Mesh& mesh);

  /** Write field over cells to file.
   *
   * @param t Time associated with field.
   * @param field Field over cells.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  void writeCellField(const PylithScalar t,
		      topology::Field& field,
		      const char* label =0,
		      const int labelId =0);

  /** Write dataset with names of points to file.
   *
   * @param names Array with name for each point, e.g., station name.
   * @param nunNames Number of names in array.
   *
   * Primarily used with OutputSolnPoints.
   */
  void writePointNames(const char* const* names,
		       const int numNames);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Copy constructor.
   *
   * @param w Object to copy.
   */
  DataWriterHDF5Ext(const DataWriterHDF5Ext& w);

  /// Generate filename for HDF5 file.
  std::string _hdf5Filename(void) const;

  /// Generate filename for external dataset file.
  std::string _datasetFilename(const char* field) const;

  /** Write time stamp to file.
   *
   * @param t Time in seconds.
   */
  void _writeTimeStamp(const PylithScalar t);  

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  const DataWriterHDF5Ext& operator=(const DataWriterHDF5Ext&); ///< Not implemented

// PRIVATE STRUCTS //////////////////////////////////////////////////////
private :

  struct ExternalDataset {
    PetscViewer viewer;
    PetscInt numTimeSteps;
    PetscInt numPoints;
    PetscInt fiberDim;
  };
  typedef std::map<std::string, ExternalDataset> dataset_type;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::string _filename; ///< Name of HDF5 file.
  HDF5* _h5; ///< HDF5 file
  dataset_type _datasets; ///< Datasets
  int _tstampIndex; ///< Index of last time stamp written.

}; // DataWriterHDF5Ext

#include "DataWriterHDF5Ext.icc" // inline methods

#endif // pylith_meshio_datawriterhdf5ext_hh


// End of file 
