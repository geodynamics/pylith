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

#include "DataWriter.hh" // USES DataWriter::FieldEnum

#include "pylith/utils/sievetypes.hh" // USES ALE::Mesh, real_section_type

namespace pylith {
  namespace meshio {
    class OutputManager;

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

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  OutputManager(void);

  /// Destructor
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
   */
  void open(const ALE::Obj<ALE::Mesh>& mesh,
	    const spatialdata::geocoords::CoordSys* csMesh,
	    const int numTimeSteps);

  /// Close output files.
  void close(void);

  /** Setup file for writing fields at time step.
   *
   * @param t Time of time step.
   * @param mesh PETSc mesh object.
   * @param csMesh Coordinate system of mesh geometry
   */
  void openTimeStep(const double t,
	       const ALE::Obj<ALE::Mesh>& mesh,
	       const spatialdata::geocoords::CoordSys* csMesh);

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
			 const DataWriter::FieldEnum fieldType,
			 const ALE::Obj<ALE::Mesh>& mesh);

  /** Append finite-element cell field to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field Cell field.
   * @param fieldType Type of field.
   * @param mesh PETSc mesh object.
   */
  void appendCellField(const double t,
		       const char* name,
		       const ALE::Obj<real_section_type>& field,
		       const DataWriter::FieldEnum fieldType,
		       const ALE::Obj<ALE::Mesh>& mesh);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  spatialdata::geocoords::CoordSys* _coordsys; ///< Coordinate system for output.
  DataWriter* _writer; ///< Writer for data.
  VertexFilter* _vertexFilter; ///< Filter applied to vertex data.
  CellFilter* _cellFilter; ///< Filter applied to cell data.

  bool _isInfo; ///< Is output info (diagnostic stuff) or data (solution, etc).

}; // OutputManager

#endif // pylith_meshio_outputmanager_hh

// End of file 
