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
 * @file pylith/meshio/OutputManager.hh
 *
 * @brief C++ object for manager output of finite-element data.
 */

#if !defined(pylith_meshio_outputmanager_hh)
#define pylith_meshio_outputmanager_hh



#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh, real_section_type
#include "pylith/utils/vectorfields.hh" // USES VectorFieldEnum

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "DataWriter.hh" // USES DataWriter in templated methods

// OutputManager --------------------------------------------------------
class pylith::meshio::OutputManager
{ // OutputManager
  friend class TestOutputManager; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  OutputManager(void);

  /// Destructor
  virtual
  ~OutputManager(void);

  /** Set coordinate system in output. The vertex fields in the output
   * are not affected by any change in coordinates.
   *
   * @param cs Coordinate system in output.
   */
  void coordsys(const spatialdata::geocoords::CoordSys* cs);

  /** Set writer to write data to file.
   *
   * @param datawriter Writer for data.
   */
  void writer(const DataWriter<mesh_type>* datawriter);

  /** Set filter for vertex data.
   *
   * @param filter Filter to apply to vertex data before writing.
   */
  void vertexFilter(const VertexFilter* filter);

  /** Set filter for cell data.
   *
   * @param filter Filter to apply to cell data before writing.
   */
  void cellFilter(const CellFilter* filter);

  /** Prepare for output.
   *
   * @param mesh Finite-element mesh object.
   * @param numTimeSteps Expected number of time steps.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  template<typename mesh_type>
  void open(const mesh_type& mesh,
	    const int numTimeSteps,
	    const char* label =0,
	    const int labelId =0);

  /// Close output files.
  void close(void);

  /** Setup file for writing fields at time step.
   *
   * @param t Time of time step.
   * @param mesh Finite-element mesh object.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  template<typename mesh_type>
  void openTimeStep(const double t,
		    const mesh_type& mesh,
		    const char* label =0,
		    const int labelId =0);

  /// End writing fields at time step.
  void closeTimeStep(void);

  /** Append finite-element vertex field to file.
   *
   * @param t Time associated with field.
   * @param field Vertex field.
   */
  template<typename mesh_type>
  void appendVertexField(const double t,
			 const topology::Field<mesh_type>& field);

  /** Append finite-element cell field to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field Cell field.
   * @param fieldType Type of field.
   * @param mesh PETSc mesh object.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  template<typename mesh_type>
  void appendCellField(const double t,
		       const char* name,
		       const topology::Field<mesh_type>& field,
		       const char* label =0,
		       const int labelId =0);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  OutputManager(const OutputManager&); ///< Not implemented.
  const OutputManager& operator=(const OutputManager&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  spatialdata::geocoords::CoordSys* _coordsys; ///< Coordinate system for output.
  DataWriter<mesh_type>* _writer; ///< Writer for data.
  VertexFilter* _vertexFilter; ///< Filter applied to vertex data.
  CellFilter* _cellFilter; ///< Filter applied to cell data.

}; // OutputManager

#include "OutputManager.icc" // template methods

#endif // pylith_meshio_outputmanager_hh


// End of file 
