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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/DataWriterVTK.hh
 *
 * @brief Object for writing finite-element data to VTK file.
 */

#if !defined(pylith_meshio_datawritervtk_hh)
#define pylith_meshio_datawritervtk_hh

// Include directives ---------------------------------------------------
#include "DataWriter.hh" // ISA DataWriter

// DataWriterVTK --------------------------------------------------------
/// Object for writing finite-element data to VTK file.
template<typename mesh_type, typename field_type>
class pylith::meshio::DataWriterVTK : public DataWriter<mesh_type,field_type>
{ // DataWriterVTK
  friend class TestDataWriterVTKMesh; // unit testing
  friend class TestDataWriterVTKSubMesh; // unit testing
  friend class TestDataWriterVTKBCMesh; // unit testing
  friend class TestDataWriterVTKFaultMesh; // unit testing
  friend class TestDataWriterVTKPoints; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  DataWriterVTK(void);

  /// Destructor
  ~DataWriterVTK(void);

  /** Make copy of this object.
   *
   * @returns Copy of this.
   */
  DataWriter<mesh_type, field_type>* clone(void) const;

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Set filename for VTK file.
   *
   * @param filename Name of VTK file.
   */
  void filename(const char* filename);

  /** Set time format for time stamp in name of VTK file.
   *
   * @param format C style time format for filename.
   */
  void timeFormat(const char* format);

  /** Set value used to normalize time stamp in name of VTK file.
   *
   * Time stamp is divided by this value (time in seconds).
   *
   * @param value Value (time in seconds) used to normalize time stamp in
   * filename.
   */
  void timeConstant(const PylithScalar value);

  /** Set precision of floating point values in output.
   *
   * @param value Precision for floating point values.
   */
  void precision(const int value);

  /** Prepare for writing files.
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

  /** Prepare file for data at a new time step.
   *
   * @param t Time stamp for new data
   * @param mesh Finite-element mesh.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  void openTimeStep(const PylithScalar t,
		    const mesh_type& mesh,
		    const char* label =0,
		    const int labelId =0);

  /// Cleanup after writing data for a time step.
  void closeTimeStep(void);

  /** Write field over vertices to file.
   *
   * @param t Time associated with field.
   * @param field Field over vertices.
   * @param mesh Mesh associated with output.
   */
  void writeVertexField(const PylithScalar t,
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
  void writeCellField(const PylithScalar t,
		      field_type& field,
		      const char* label =0,
		      const int labelId =0);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Copy constructor.
   *
   * @param w Object to copy.
   */
  DataWriterVTK(const DataWriterVTK& w);

  /** Generate filename for VTK file.
   *
   * @param t Time in seconds.
   */
  std::string _vtkFilename(const PylithScalar t) const;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  const DataWriterVTK& operator=(const DataWriterVTK&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /// Time value (in seconds) used to normalize time stamp.
  PylithScalar _timeConstant;

  std::string _filename; ///< Name of VTK file.
  std::string _timeFormat; ///< C style time format for time stamp.

  PetscViewer _viewer; ///< Output file
  PetscDM _dm; ///< Handle to PETSc DM for mesh

  int _precision; ///< Precision of floating point values in output.

  bool _isOpen; ///< True if called open().
  bool _isOpenTimeStep; ///< true if called openTimeStep().
  bool _wroteVertexHeader; ///< True if wrote header for vertex data.
  bool _wroteCellHeader; ///< True if wrote header for cell data

}; // DataWriterVTK

#include "DataWriterVTK.icc" // inline methods
#include "DataWriterVTK.cc" // template definitions

#endif // pylith_meshio_datawritervtk_hh


// End of file 
