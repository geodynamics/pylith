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

#include "pylith/utils/sievetypes.hh" // USES ALE::Obj, ALE::Mesh, real_section_type

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
   */
  virtual
  void open(const ALE::Obj<ALE::Mesh>& mesh,
	    const spatialdata::geocoords::CoordSys* csMesh);

  /// Close output files.
  virtual
  void close(void);

  /** Prepare file for data at a new time step.
   *
   * @param t Time stamp for new data
   * @param mesh PETSc mesh object
   * @param csMesh Coordinate system of mesh geometry
   */
  virtual
  void openTimeStep(const double t,
		    const ALE::Obj<ALE::Mesh>& mesh,
		    const spatialdata::geocoords::CoordSys* csMesh);

  /// Cleanup after writing data for a time step.
  virtual
  void closeTimeStep(void);

  /** Write field over vertices to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field PETSc field over vertices.
   * @param mesh Finite-element mesh
   */
  virtual
  void writeVertexField(const double t,
			const char* name,
			const ALE::Obj<real_section_type>& field,
			const ALE::Obj<ALE::Mesh>& mesh) = 0;

  /** Write field over cells to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field PETSc field over cells.
   * @param mesh PETSc mesh object.
   */
  virtual
  void writeCellField(const double t,
		      const char* name,
		      const ALE::Obj<real_section_type>& field,
		      const ALE::Obj<ALE::Mesh>& mesh) = 0;

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

}; // DataWriter

#endif // pylith_meshio_datawriter_hh


// End of file 
