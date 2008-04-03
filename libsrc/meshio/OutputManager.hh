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

namespace pylith {
  namespace meshio {
    class OutputManager;
    class TestOutputManager; // unit testing

    class DataWriter; // HOLDS DataWriter
    class CellFilter; // HOLDSA CellFilter
    class VertexFilter; // HOLDSA VertexFilter
  } // meshio
} // pylith

namespace spatialdata {
  namespace geocoords {
    class CoordSys; // USES CoordSys
  } // geocoords
} // spatialdata

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
  void writer(const DataWriter* datawriter);

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
   * @param mesh PETSc mesh object.
   * @param csMesh Coordinate system of mesh geometry.
   * @param numTimeSteps Expected number of time steps.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  void open(const ALE::Obj<Mesh>& mesh,
	    const spatialdata::geocoords::CoordSys* csMesh,
	    const int numTimeSteps,
	    const char* label =0,
	    const int labelId =0);

  /// Close output files.
  void close(void);

  /** Setup file for writing fields at time step.
   *
   * @param t Time of time step.
   * @param mesh PETSc mesh object.
   * @param csMesh Coordinate system of mesh geometry
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  void openTimeStep(const double t,
		    const ALE::Obj<Mesh>& mesh,
		    const spatialdata::geocoords::CoordSys* csMesh,
		    const char* label =0,
		    const int labelId =0);

  /// End writing fields at time step.
  void closeTimeStep(void);

  /** Append finite-element vertex field to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field Vertex field.
   * @param fieldType Type of field.
   * @param mesh PETSc mesh object.
   */
  void appendVertexField(const double t,
			 const char* name,
			 const ALE::Obj<real_section_type>& field,
			 const VectorFieldEnum fieldType,
			 const ALE::Obj<Mesh>& mesh);

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
  void appendCellField(const double t,
		       const char* name,
		       const ALE::Obj<real_section_type>& field,
		       const VectorFieldEnum fieldType,
		       const ALE::Obj<Mesh>& mesh,
		       const char* label =0,
		       const int labelId =0);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  OutputManager(const OutputManager&); ///< Not implemented.
  const OutputManager& operator=(const OutputManager&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  spatialdata::geocoords::CoordSys* _coordsys; ///< Coordinate system for output.
  DataWriter* _writer; ///< Writer for data.
  VertexFilter* _vertexFilter; ///< Filter applied to vertex data.
  CellFilter* _cellFilter; ///< Filter applied to cell data.

}; // OutputManager

#endif // pylith_meshio_outputmanager_hh

// End of file 
