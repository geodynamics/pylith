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

#include "pylith/utils/sievetypes.hh" // USES ALE::Obj, ALE::Mesh, real_section_type

namespace pylith {
  namespace meshio {
    class OutputManager;

    class DataWriter; // HOLDS DataWriter
    class VertexFilter; // HOLDSA VertexFilter
    class CellFilter; // HOLDSA CellFilter
  } // meshio

  namespace topology {
    class FieldsManager; // USES FieldsManager
  } // topology
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

  /** Set writer to write data to file.
   *
   * @param datawriter Writer for data.
   */
  void writer(const DataWriter* datawriter);

  /** Set which vertex fields to output.
   *
   * @param names Names of fields.
   * @param nfields Number of fields
   */
  void vertexFields(const char** names,
		    const int nfields);

  /** Set which cell fields to output.
   *
   * @param names Names of fields.
   * @param nfields Number of fields
   */
  void cellFields(const char** names,
		  const int nfields);

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
   * @param mesh PETSc mesh object
   * @param csMesh Coordinate system of mesh geometry
   */
  void open(const ALE::Obj<ALE::Mesh>& mesh,
	    const spatialdata::geocoords::CoordSys* csMesh);

  /// Close output files.
  void close(void);

  /** Write finite-element fields to file.
   *
   * @param t Time associated with field.
   * @param fields Fields manager.
   * @param mesh PETSc mesh object.
   * @param csMesh Coordinate system of mesh geometry
   */
  void writeFields(const double t,
		   const topology::FieldsManager* fields,
		   const ALE::Obj<ALE::Mesh>& mesh,
		   const spatialdata::geocoords::CoordSys* csMesh);

// PRIVATE MEMBERS //////////////////////////////////////////////////////

private :

  std::vector<std::string> _vertexFields; ///< Names of vertex fields to output
  std::vector<std::string> _cellFields; ///< Names of cell fields to output

  DataWriter* _writer; ///< Writer for data
  VertexFilter* _vertexFilter; ///< Filter applied to vertex data
  CellFilter* _cellFilter; ///< Filter applied to cell data

}; // OutputManager

#endif // pylith_meshio_outputmanager_hh

// End of file 
