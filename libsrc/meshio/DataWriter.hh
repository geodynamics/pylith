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

/**
 * @file pylith/meshio/DataWriter.hh
 *
 * @brief Abstract base class for writing finite-element data to file.
 */

#if !defined(pylith_meshio_datawriter_hh)
#define pylith_meshio_datawriter_hh

#include "pylith/utils/sievetypes.hh" // USES ALE::Obj, PETSc Mesh, real_section_type
#include "pylith/utils/vectorfields.hh" // USES VectorFieldEnum

namespace pylith {
  namespace meshio {
    class DataWriter;
  } // meshio
} // pylith

namespace spatialdata {
  namespace geocoords {
    class CoordSys; // USES CoordSys
  } // geocoords
} // spatialdata

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

  /** Prepare for writing files.
   *
   * @param mesh PETSc mesh object 
   * @param csMesh Coordinate system of mesh geometry
   * @param numTimeSteps Expected number of time steps for fields.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  virtual
  void open(const ALE::Obj<Mesh>& mesh,
	    const spatialdata::geocoords::CoordSys* csMesh,
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
   * @param csMesh Coordinate system of mesh geometry
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  virtual
  void openTimeStep(const double t,
		    const ALE::Obj<Mesh>& mesh,
		    const spatialdata::geocoords::CoordSys* csMesh,
		    const char* label =0,
		    const int labelId =0);

  /// Cleanup after writing data for a time step.
  virtual
  void closeTimeStep(void);

  /** Write field over vertices to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field PETSc field over vertices.
   * @param fieldType Type of field.
   * @param mesh Finite-element mesh
   */
  virtual
  void writeVertexField(const double t,
			const char* name,
			const ALE::Obj<real_section_type>& field,
			const VectorFieldEnum fieldType,
			const ALE::Obj<Mesh>& mesh) = 0;

  /** Write field over cells to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field PETSc field over cells.
   * @param fieldType Type of field.
   * @param mesh PETSc mesh object.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  virtual
  void writeCellField(const double t,
		      const char* name,
		      const ALE::Obj<real_section_type>& field,
		      const VectorFieldEnum fieldType,
		      const ALE::Obj<Mesh>& mesh,
		      const char* label =0,
		      const int labelId =0) = 0;

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

  int _numTimeSteps; ///< Expected number of time steps for fields.

}; // DataWriter

#endif // pylith_meshio_datawriter_hh


// End of file 
