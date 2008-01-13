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

#include "pylith/utils/sievetypes.hh" // USES ALE::Mesh, real_section_type
#include <map> // USES std::map

namespace pylith {
  namespace meshio {
    class OutputManager;

    class DataWriter; // HOLDS DataWriter
    class CellFilter; // HOLDSA CellFilter
    class VertexFilter; // HOLDSA VertexFilter
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

  /// Map to hold field names and mesh labels (name -> label).
  typedef std::map<std::string, std::string> map_names_type;

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
   * @param labels Mesh labels of fields.
   * @param numFields Number of fields.
   */
  void vertexFields(const char** names,
		    const char** labels,
		    const int numFields);

  /** Set which cell fields to output.
   *
   * @param names Names of fields.
   * @param labels Mesh labels of fields.
   * @param numFields Number of fields.
   */
  void cellFields(const char** names,
		  const char** labels,
		  const int numFields);

  /** Get vertex fields to output.
   *
   * @returns Map of field name to mesh label for fields.
   */
  const map_names_type& vertexFields(void) const;

  /** Get cell fields to output.
   *
   * @returns Map of field name to mesh label for fields.
   */
  const map_names_type& cellFields(void) const;

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
		   topology::FieldsManager* const fields,
		   const ALE::Obj<ALE::Mesh>& mesh,
		   const spatialdata::geocoords::CoordSys* csMesh);

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
   * @param mesh PETSc mesh object.
   * @param csMesh Coordinate system of mesh geometry
   */
  void appendVertexField(const double t,
			 const char* name,
			 const ALE::Obj<real_section_type>& field,
			 const ALE::Obj<ALE::Mesh>& mesh,
			 const spatialdata::geocoords::CoordSys* csMesh);

  /** Append finite-element cell field to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field Cell field.
   * @param mesh PETSc mesh object.
   * @param csMesh Coordinate system of mesh geometry
   */
  void appendCellField(const double t,
		       const char* name,
		       const ALE::Obj<real_section_type>& field,
		       const ALE::Obj<ALE::Mesh>& mesh,
		       const spatialdata::geocoords::CoordSys* csMesh);

// PRIVATE MEMBERS //////////////////////////////////////////////////////

private :

  /// Name and section label of vertex fields to output
  map_names_type _vertexFields;

  /// Name and section label of cell fields to output
  map_names_type _cellFields;

  DataWriter* _writer; ///< Writer for data
  VertexFilter* _vertexFilter; ///< Filter applied to vertex data
  CellFilter* _cellFilter; ///< Filter applied to cell data

}; // OutputManager

#endif // pylith_meshio_outputmanager_hh

// End of file 
