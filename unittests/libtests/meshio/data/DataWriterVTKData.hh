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

#if !defined(pylith_meshio_datawritervtkdata_hh)
#define pylith_meshio_datawritervtkdata_hh

#include "pylith/topology/FieldBase.hh" // USES VectorFieldEnum

namespace pylith {
  namespace meshio {
     class DataWriterVTKData;
  } // meshio
} // pylith

class pylith::meshio::DataWriterVTKData
{ // DataWriterVTKData

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  DataWriterVTKData(void);

  /// Destructor
  virtual
  ~DataWriterVTKData(void);

// PUBLIC STRUCTS ///////////////////////////////////////////////////////
public:

  struct FieldStruct {
    char* name; ///< Name of field
    topology::FieldBase::VectorFieldEnum field_type; ///< Type of field.
    int fiber_dim; ///< Fiber dimension for field.
  }; // FieldStruct

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* meshFilename; ///< Name of mesh file.
  char* faultLabel; ///< Name of group of vertices for fault.
  int faultId; ///< Material identifier for fault.
  char* bcLabel; ///< Name of group of vertices for bc.

  char* timestepFilename; ///< Name of VTK file without fields.
  char* vertexFilename; ///< Name of VTK file for vertex fields.
  char* cellFilename; ///< Name of VTK file for cell fields.

  double time; ///< Time for fields.
  char* timeFormat; ///< Format for time stamp.

  char* cellsLabel; ///< Name of label for mesh cells (if using subset or boundary).
  int labelId; ///< Id for label associated with cells (if cellsLabel != 0)

  /// @name Vertex field information.
  //@{
  int numVertexFields; ///< Number of vertex fields.
  int numVertices; ///< Number of vertices.
  FieldStruct* vertexFieldsInfo; ///< Array of vertex field information.
  double* vertexFields[3]; ///< Array of vertex field values.
  //@}

  /// @name Cell field information.
  //@{
  int numCellFields; ///< Number of cell fields.
  int numCells; ///< Number of vertices.
  FieldStruct* cellFieldsInfo; ///< Array of cell fields information.
  double* cellFields[3]; /// Array of cell field values.
  //@}

}; // DataWriterVTKData

#endif // pylith_meshio_datawritervtkdata_hh


// End of file
