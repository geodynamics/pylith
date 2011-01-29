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

/**
 * @file libsrc/meshio/DataWriterHDF5.hh
 *
 * @brief Object for writing finite-element data to HDF5 file.
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

#if !defined(pylith_meshio_datawriterhdf5_hh)
#define pylith_meshio_datawriterhdf5_hh

// Include directives ---------------------------------------------------
#include "DataWriter.hh" // ISA DataWriter

#include <string> // USES std::string
#include <map> // HASA std::map

// DataWriterHDF5 --------------------------------------------------------
/// Object for writing finite-element data to HDF5 file.
template<typename mesh_type, typename field_type>
class pylith::meshio::DataWriterHDF5 : public DataWriter<mesh_type,field_type>
{ // DataWriterHDF5
  friend class TestDataWriterHDF5Mesh; // unit testing
  friend class TestDataWriterHDF5SubMesh; // unit testing
  friend class TestDataWriterHDF5BCMesh; // unit testing
  friend class TestDataWriterHDF5FaultMesh; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  DataWriterHDF5(void);

  /// Destructor
  ~DataWriterHDF5(void);

  /** Make copy of this object.
   *
   * @returns Copy of this.
   */
  DataWriter<mesh_type, field_type>* clone(void) const;

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Set filename for HDF5 file.
   *
   * @param filename Name of HDF5 file.
   */
  void filename(const char* filename);

  /** Open output file.
   *
   * @param mesh Finite-element mesh. 
   * @param numTimeSteps Expected number of time steps for fields.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  void open(const mesh_type& mesh,
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
  void writeVertexField(const double t,
			field_type& field,
			const mesh_type& mesh);

  /** Write field over cells to file.
   *
   * @param t Time associated with field.
   * @param field Field over cells.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  void writeCellField(const double t,
		      field_type& field,
		      const char* label =0,
		      const int labelId =0);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Copy constructor.
   *
   * @param w Object to copy.
   */
  DataWriterHDF5(const DataWriterHDF5& w);

  /// Generate filename for HDF5 file.
  std::string _hdf5Filename(void) const;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  const DataWriterHDF5& operator=(const DataWriterHDF5&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::string _filename; ///< Name of HDF5 file.
  PetscViewer _viewer; ///< Output file

  std::map<std::string, int> _timesteps; ///< # of time steps written per field.

}; // DataWriterHDF5

#include "DataWriterHDF5.icc" // inline methods
#include "DataWriterHDF5.cc" // template definitions

#endif // pylith_meshio_datawriterhdf5_hh


// End of file 
