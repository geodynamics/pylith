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

/**
 * @file libsrc/meshio/DataWriter.hh
 *
 * @brief Abstract base class for writing finite-element data to file.
 */

#if !defined(pylith_meshio_datawriter_hh)
#define pylith_meshio_datawriter_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Field
#include <string> // HASA std::string

// DataWriter -----------------------------------------------------------
/// Abstract base class for writing finite-element data to file.
class pylith::meshio::DataWriter
{ // DataWriter

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  DataWriter(void);

  /// Destructor
  virtual
  ~DataWriter(void);

  /** Make copy of this object.
   *
   * @returns Copy of this.
   */
  virtual
  DataWriter* clone(void) const = 0;

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set time scale for simulation time.
   *
   * @param value Time scale
   */
  void timeScale(const PylithScalar value);

  /** Prepare for writing files.
   *
   * @param mesh Finite-element mesh. 
   * @param numTimeSteps Expected number of time steps for fields.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  virtual
  void open(const topology::Mesh& mesh,
	    const int numTimeSteps,
	    const char* label =0,
	    const int labelId =0);

  /// Close output files.
  virtual
  void close(void);

  /** Prepare file for data at a new time step.
   *
   * @param t Time stamp for new data
   * @param mesh PETSc mesh object
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  virtual
  void openTimeStep(const PylithScalar t,
		    const topology::Mesh& mesh,
		    const char* label =0,
		    const int labelId =0);

  /// Cleanup after writing data for a time step.
  virtual
  void closeTimeStep(void);

  /** Write field over vertices to file.
   *
   * @param t Time associated with field.
   * @param field Field over vertices.
   * @param mesh Mesh associated with output.
  */
  virtual
  void writeVertexField(const PylithScalar t,
			topology::Field& field,
			const topology::Mesh& mesh) = 0;

  /** Write field over cells to file.
   *
   * @param t Time associated with field.
   * @param field Field over cells.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  virtual
  void writeCellField(const PylithScalar t,
		      topology::Field& field,
		      const char* label =0,
		      const int labelId =0) = 0;

  /** Write dataset with names of points to file.
   *
   * @param names Array with name for each point, e.g., station name.
   * @param nunNames Number of names in array.
   * @param mesh Finite-element mesh. 
   *
   * Primarily used with OutputSolnPoints.
   */
  virtual
  void writePointNames(const char* const* names,
		       const int numNames,
		       const topology::Mesh& mesh);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param w Object to copy.
   */
  DataWriter(const DataWriter& w);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  const DataWriter& operator=(const DataWriter&); ///< Not implemented

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  PylithScalar _timeScale; ///< Time scale for dimensioning time in output.
  int _numTimeSteps; ///< Expected number of time steps for fields.
  std::string _context; ///< Context of scatters for DataWriter.

}; // DataWriter

#endif // pylith_meshio_datawriter_hh


// End of file 
